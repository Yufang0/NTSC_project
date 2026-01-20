############################################################
## 0) 套件
############################################################
pkgs <- c(
  "survival", "dplyr", "Matrix",
  "glmnet",         # LASSO Cox
  "grpreg",         # Group Lasso (selection step)
  "randomForestSRC",# RSF
  "xgboost",        # XGBoost survival:cox
  "SHAPforxgboost"  # SHAP for xgboost
)
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

############################################################
## 1) 載入資料（預設用 survival 套件內建 gbsg）
############################################################
LOAD_DATA <- "pkg"  # "pkg" or "csv"

if (LOAD_DATA == "pkg") {
  data(gbsg, package = "survival")  # 內建資料
  df <- gbsg
  # survival::gbsg 欄位通常是：
  # age, meno(0/1), size, grade, nodes, pgr, er, hormon(0/1), rfstime, status(0/1), pid
} else {
  # 若你要用你上傳的 CSV，請先重新上傳檔案並改路徑
  df <- read.csv("/mnt/data/gbsg.csv", stringsAsFactors = FALSE)
  
  # 若出現 X / X.1 之類索引欄，移除
  if ("X" %in% names(df)) df <- df %>% dplyr::select(-X)
  if ("X.1" %in% names(df)) df <- df %>% dplyr::select(-`X.1`)
}

############################################################
## 2) 基本清理/型別
############################################################
df <- df %>%
  mutate(
    # meno/hormon 可能是 0/1 或 pre/post
    meno   = if (is.numeric(meno)) factor(ifelse(meno == 1, "post", "pre")) else factor(meno),
    hormon = if (is.numeric(hormon)) factor(ifelse(hormon == 1, "yes", "no")) else factor(hormon),
    grade  = factor(grade),
    status = as.integer(status),
    rfstime = as.numeric(rfstime)
  )

# 移除 pid（識別碼不當特徵）
if ("pid" %in% names(df)) df <- df %>% dplyr::select(-pid)

# 去 NA
df <- na.omit(df)

############################################################
## 3) Outer Train/Test split 70/30 + Inner split（避免洩漏）
############################################################
set.seed(2025)
n <- nrow(df)

# ---- outer split：最後一次性評估用（test_outer 永遠不拿來 fit / 選參數）----
idx_outer_tr <- sample(seq_len(n), size = floor(0.7 * n))
train_outer <- df[idx_outer_tr, ]
test_outer  <- df[-idx_outer_tr, ]

# ---- inner split：所有「調參/選變數/比較」都在 train_outer 裡面完成----
set.seed(2025)
n_tr <- nrow(train_outer)
idx_inner_tr <- sample(seq_len(n_tr), size = floor(0.8 * n_tr))  # 80/20 你也可改 70/30

train_inner <- train_outer[idx_inner_tr, ]
val_inner   <- train_outer[-idx_inner_tr, ]

# 建立 Surv 物件（後面用）
y_inner <- Surv(train_inner$rfstime, train_inner$status)
y_val   <- Surv(val_inner$rfstime,   val_inner$status)
y_test  <- Surv(test_outer$rfstime,  test_outer$status)

############################################################
## 4) 評估函數：C-index + (事件者) MSE(預測中位存活時間)
##    注意：MSE 對 censoring 不是最嚴謹，但你流程圖寫了 MSE，
##          這裡用「只對 status=1 的人」做示範。
############################################################
c_index <- function(time, status, risk_score) {
  # risk_score 越大 = 越高風險（越早事件）
  # concordance 預設方向是「分數越大越長壽」，所以要反向
  survival::concordance(Surv(time, status) ~ I(-risk_score))$concordance
}


median_from_survcurve <- function(time_grid, surv_prob) {
  # 找到 S(t) <= 0.5 的最早時間
  k <- which(surv_prob <= 0.5)
  if (length(k) == 0) return(NA_real_)
  time_grid[min(k)]
}

mse_median_time_on_events <- function(time, status, pred_median_time) {
  ev <- which(status == 1L & !is.na(pred_median_time))
  if (length(ev) == 0) return(NA_real_)
  mean((pred_median_time[ev] - time[ev])^2)
}

############################################################
## 5) (A) 隨機森林：Random Survival Forest (RSF)
############################################################
set.seed(2025)
rsf_fit <- randomForestSRC::rfsrc(
  Surv(rfstime, status) ~ .,
  data = train_outer,
  ntree = 1500,
  nodesize = 15,
  importance = "permute"
)

# RSF 風險分數（越大通常越高風險；不同版本欄位名稱略不同）
rsf_pred <- predict(rsf_fit, newdata = test_outer)

# 多數情況可用 rsf_pred$predicted 當 risk score
rsf_risk <- as.numeric(rsf_pred$predicted)

rsf_cidx <- c_index(test_outer$rfstime, test_outer$status, rsf_risk)

# RSF 若要拿 survival curve 做中位數（版本差異大），這裡示範用風險排序即可
rsf_mse  <- NA_real_

############################################################
## 6) (B) XGBoost survival:cox + SHAP
############################################################
# one-hot encoding
x_tr <- model.matrix(~ . - 1, data = train_outer %>% dplyr::select(-rfstime, -status))
x_te <- model.matrix(~ . - 1, data = test_outer  %>% dplyr::select(-rfstime, -status))

# 對齊欄位
missing_in_te <- setdiff(colnames(x_tr), colnames(x_te))
if (length(missing_in_te) > 0) {
  x_te <- cbind(x_te, matrix(0, nrow = nrow(x_te), ncol = length(missing_in_te),
                             dimnames = list(NULL, missing_in_te)))
}
x_te <- x_te[, colnames(x_tr), drop = FALSE]

# XGBoost survival:cox 標籤：事件用 +time，censor 用 -time :contentReference[oaicite:4]{index=4}
label_tr <- ifelse(train_outer$status == 1L, train_outer$rfstime, -train_outer$rfstime)
label_te <- ifelse(test_outer$status  == 1L, test_outer$rfstime,  -test_outer$rfstime)

dtr <- xgboost::xgb.DMatrix(data = x_tr, label = label_tr)
dte <- xgboost::xgb.DMatrix(data = x_te, label = label_te)

params <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik",
  eta = 0.03,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1
)

set.seed(2025)
xgb_fit <- xgboost::xgb.train(
  params = params,
  data = dtr,
  nrounds = 1500,
  verbose = 0
)

xgb_risk <- predict(xgb_fit, dte)  # 這是 Cox risk score（線性 predictor）
xgb_cidx <- c_index(test_outer$rfstime, test_outer$status, xgb_risk)

# SHAP（針對 xgb）
# 計算 SHAP values（用 x_tr / x_te 都可；建議用 train_outer 看特徵影響）
sh <- SHAPforxgboost::shap.values(xgb_model = xgb_fit, X_train = x_tr)
# shap_long <- SHAPforxgboost::shap.prep(shap_contrib = sh$shap_score, X_train = x_tr)
# 如要畫圖：SHAPforxgboost::shap.plot.summary(shap_long)

xgb_mse <- NA_real_

############################################################
## 7) (C) Cox (naive) + SIMSELEX（雙層 split：inner 做比較，outer test 做最終）
############################################################

# ---- helper：建立與 test_outer 相同欄位的 design matrix ----
make_xmm <- function(dat, time_col = "rfstime", status_col = "status") {
  x <- model.matrix(
    reformulate(setdiff(names(dat), c(time_col, status_col))),
    data = dat
  )
  if ("(Intercept)" %in% colnames(x)) x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  x
}

# ---- helper：用「訓練集」的 baseline cumulative hazard + lp(newdata) 算 survival curve ----
# S(t|lp) = exp( -H0(t) * exp(lp) )
surv_from_basehaz <- function(basehaz_df, lp_vec) {
  tgrid <- basehaz_df$time
  H0    <- basehaz_df$hazard
  S_mat <- exp(- outer(H0, exp(lp_vec)))
  list(time = tgrid, surv = S_mat)
}

# ---------- 7.1 Cox naive：inner 評估 ----------
cox_naive_inner <- coxph(Surv(rfstime, status) ~ ., data = train_inner, x = TRUE)

lp_val_naive <- predict(cox_naive_inner, newdata = val_inner, type = "lp")
cidx_val_naive <- c_index(val_inner$rfstime, val_inner$status, lp_val_naive)

sf_val_naive <- survfit(cox_naive_inner, newdata = val_inner)
tgrid_naive  <- sf_val_naive$time
med_val_naive <- apply(sf_val_naive$surv, 2, function(s) median_from_survcurve(tgrid_naive, s))
mse_val_naive <- mse_median_time_on_events(val_inner$rfstime, val_inner$status, med_val_naive)

# ---------- 7.2 Cox naive：outer 最終評估（在 train_outer 重跑一次） ----------
cox_naive_outer <- coxph(Surv(rfstime, status) ~ ., data = train_outer, x = TRUE)

lp_test_naive <- predict(cox_naive_outer, newdata = test_outer, type = "lp")
cox_cidx <- c_index(test_outer$rfstime, test_outer$status, lp_test_naive)

sf_test_naive <- survfit(cox_naive_outer, newdata = test_outer)
tgrid_test_naive <- sf_test_naive$time
cox_med <- apply(sf_test_naive$surv, 2, function(s) median_from_survcurve(tgrid_test_naive, s))
cox_mse <- mse_median_time_on_events(test_outer$rfstime, test_outer$status, cox_med)

# ---------- 7.3 SIMSELEX：函數（沿用你原本，僅保留 sqrt(lambda) + group-lasso + extrapolate） ----------
simselex_cox <- function(dat, time_col, status_col,
                         kappa = 0.8,
                         error_vars = c("size"),
                         lambdas = c(0, 0.5, 1, 1.5, 2),
                         B = 50,
                         poly_degree = 2,
                         lasso_s = "lambda.min",
                         clamp_nonneg = TRUE,
                         seed = 2025) {
  
  set.seed(seed)
  
  y <- Surv(dat[[time_col]], dat[[status_col]])
  x_mm <- make_xmm(dat, time_col, status_col)
  
  var_names <- colnames(x_mm)
  err_idx <- which(var_names %in% error_vars)
  if (length(err_idx) == 0) stop("error_vars 沒有在 model.matrix 後欄名中找到。請檢查名稱。")
  
  # kappa = var(X)/(var(X)+var(d)) -> var(d)=var(X)*(1-kappa)/kappa
  sigma_d2 <- apply(x_mm[, err_idx, drop = FALSE], 2, var) * (1 - kappa) / kappa
  sigma_d  <- sqrt(pmax(sigma_d2, 0))
  
  K <- length(lambdas)
  p <- ncol(x_mm)
  n <- nrow(x_mm)
  n_err <- length(err_idx)
  
  beta_bar <- matrix(0, nrow = K, ncol = p,
                     dimnames = list(paste0("lam_", lambdas), var_names))
  
  for (k in seq_len(K)) {
    lam <- lambdas[k]
    beta_b <- matrix(0, nrow = B, ncol = p)
    
    for (b in seq_len(B)) {
      x_star <- x_mm
      
      if (lam > 0) {
        noise <- matrix(rnorm(n * n_err), nrow = n, ncol = n_err)
        noise <- sweep(noise, 2, sigma_d, `*`)             # 每欄自己的 sd
        x_star[, err_idx] <- x_star[, err_idx] + sqrt(lam) * noise
        
        if (clamp_nonneg) {
          x_star[, err_idx] <- pmax(x_star[, err_idx], 0)
        }
      }
      
      cvfit <- glmnet::cv.glmnet(
        x = x_star, y = y,
        family = "cox", alpha = 1,
        nfolds = 5
      )
      beta_b[b, ] <- as.numeric(as.matrix(coef(cvfit, s = lasso_s)))
    }
    
    beta_bar[k, ] <- colMeans(beta_b)
  }
  
  # selection：group lasso on polynomial coefficients
  q <- poly_degree + 1
  Z <- sapply(0:poly_degree, function(d) lambdas^d)  # K x q
  
  y_vec <- as.numeric(t(beta_bar))  # length p*K
  
  X_big <- matrix(0, nrow = p * K, ncol = p * q)
  grp   <- rep(seq_len(p), each = q)
  
  row_id <- 1
  for (j in seq_len(p)) {
    cols_j <- ((j - 1) * q + 1):((j - 1) * q + q)
    for (kk in seq_len(K)) {
      X_big[row_id, cols_j] <- Z[kk, ]
      row_id <- row_id + 1
    }
  }
  
  set.seed(seed)
  cvg <- grpreg::cv.grpreg(X_big, y_vec, group = grp, penalty = "grLasso", family = "gaussian")
  
  gamma_hat <- as.numeric(coef(cvg$fit, lambda = cvg$lambda.min))[-1]  # drop intercept
  gamma_mat <- matrix(gamma_hat, nrow = q, ncol = p)
  
  selected <- which(colSums(abs(gamma_mat)) > 0)
  
  # extrapolate to lambda = -1
  z_minus1 <- (-1)^(0:poly_degree)
  beta_simselex <- rep(0, p)
  for (j in selected) beta_simselex[j] <- sum(z_minus1 * gamma_mat[, j])
  names(beta_simselex) <- var_names
  
  list(
    x_mm = x_mm,
    beta_bar = beta_bar,
    selected = var_names[selected],
    beta_simselex = beta_simselex
  )
}

error_vars <- c("size")  # 你可擴到 c("size","nodes","pgr","er")

# ---------- 7.4 SIMSELEX：inner 評估 ----------
simsel_inner <- simselex_cox(
  dat = train_inner,
  time_col = "rfstime",
  status_col = "status",
  kappa = 0.8,
  error_vars = error_vars,
  lambdas = c(0, 0.5, 1, 1.5, 2),
  B = 50,
  poly_degree = 2,
  clamp_nonneg = TRUE,
  lasso_s = "lambda.min"
)

x_val_mm <- make_xmm(val_inner, "rfstime", "status")
# 對齊欄位（以 train_inner 的欄位為準）
missing_val <- setdiff(colnames(simsel_inner$x_mm), colnames(x_val_mm))
if (length(missing_val) > 0) {
  x_val_mm <- cbind(x_val_mm,
                    matrix(0, nrow = nrow(x_val_mm), ncol = length(missing_val),
                           dimnames = list(NULL, missing_val)))
}
x_val_mm <- x_val_mm[, colnames(simsel_inner$x_mm), drop = FALSE]

beta_sim_inner <- simsel_inner$beta_simselex
lp_val_sim <- as.numeric(x_val_mm %*% beta_sim_inner)
cidx_val_sim <- c_index(val_inner$rfstime, val_inner$status, lp_val_sim)

# baseline hazard 只能用 train_inner 估
lp_train_inner <- as.numeric(simsel_inner$x_mm %*% beta_sim_inner)
cox_offset_inner <- coxph(Surv(rfstime, status) ~ offset(lp_train_inner), data = train_inner, x = TRUE)
bh_inner <- basehaz(cox_offset_inner, centered = FALSE)

sc_val <- surv_from_basehaz(bh_inner, lp_val_sim)
med_val_sim <- apply(sc_val$surv, 2, function(s) median_from_survcurve(sc_val$time, s))
mse_val_sim <- mse_median_time_on_events(val_inner$rfstime, val_inner$status, med_val_sim)

# ---------- 7.5 SIMSELEX：outer 最終評估（在 train_outer 重跑一次） ----------
simsel_outer <- simselex_cox(
  dat = train_outer,
  time_col = "rfstime",
  status_col = "status",
  kappa = 0.8,
  error_vars = error_vars,
  lambdas = c(0, 0.5, 1, 1.5, 2),
  B = 50,
  poly_degree = 2,
  clamp_nonneg = TRUE,
  lasso_s = "lambda.min"
)

x_test_mm <- make_xmm(test_outer, "rfstime", "status")
missing_test <- setdiff(colnames(simsel_outer$x_mm), colnames(x_test_mm))
if (length(missing_test) > 0) {
  x_test_mm <- cbind(x_test_mm,
                     matrix(0, nrow = nrow(x_test_mm), ncol = length(missing_test),
                            dimnames = list(NULL, missing_test)))
}
x_test_mm <- x_test_mm[, colnames(simsel_outer$x_mm), drop = FALSE]

beta_sim_outer <- simsel_outer$beta_simselex
lp_test_sim <- as.numeric(x_test_mm %*% beta_sim_outer)
sim_cidx <- c_index(test_outer$rfstime, test_outer$status, lp_test_sim)

# baseline hazard 只能用 train_outer 估（嚴禁用 test_outer）
lp_train_outer <- as.numeric(simsel_outer$x_mm %*% beta_sim_outer)
cox_offset_outer <- coxph(Surv(rfstime, status) ~ offset(lp_train_outer), data = train_outer, x = TRUE)
bh_outer <- basehaz(cox_offset_outer, centered = FALSE)

sc_test <- surv_from_basehaz(bh_outer, lp_test_sim)
sim_med <- apply(sc_test$surv, 2, function(s) median_from_survcurve(sc_test$time, s))
sim_mse <- mse_median_time_on_events(test_outer$rfstime, test_outer$status, sim_med)

# 你可以把 inner 的結果印出來，做為「內部比較/選參數」的依據
cat("\n[Inner validation] Cox naive: C-index=", cidx_val_naive, " MSE=", mse_val_naive, "\n")
cat("[Inner validation] Cox SIMSELEX: C-index=", cidx_val_sim, " MSE=", mse_val_sim, "\n")

cat("\n[Outer test] Cox naive: C-index=", cox_cidx, " MSE=", cox_mse, "\n")
cat("[Outer test] Cox SIMSELEX: C-index=", sim_cidx, " MSE=", sim_mse, "\n")

cat("\nSIMSELEX selected variables (outer fit):\n")
print(simsel_outer$selected)

############################################################
## 8) 統整輸出（照你流程圖：C-index + MSE）
############################################################
result <- data.frame(
  Model = c("RSF", "XGBoost(survival:cox)", "Cox(naive)", "Cox(SIMSELEX)"),
  C_index = c(rsf_cidx, xgb_cidx, cox_cidx, sim_cidx),
  MSE_median_time_events_only = c(rsf_mse, xgb_mse, cox_mse, sim_mse)
)

print(result)
cat("\nSIMSELEX selected variables (outer fit):\n")
print(simsel_outer$selected)

cat("\nSIMSELEX selected variables (inner fit):\n")
print(simsel_inner$selected)

