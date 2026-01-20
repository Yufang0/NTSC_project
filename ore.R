setwd("C:/Users/Jasmin/大專生研究計畫/data/")

## =========================================================
## 0) Packages
## =========================================================
pkgs <- c("survival", "randomForestSRC", "riskRegression", "prodlim", "irr", "dplyr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(survival)
library(randomForestSRC)
library(riskRegression)
library(prodlim)
library(irr)
library(dplyr)

set.seed(2026)

## =========================================================
## 1) Read your dataset (uploaded gbsg.csv)
## =========================================================
df <- read.csv("gbsg.csv", stringsAsFactors = FALSE)

## 欄位型別（依你給的欄位定義）
df <- df %>%
  mutate(
    pid    = as.integer(pid),
    meno   = factor(meno, levels = c(0, 1), labels = c("pre", "post")),
    hormon = factor(hormon, levels = c(0, 1), labels = c("no", "yes")),
    grade  = factor(grade),
    status = as.integer(status),
    rfstime = as.numeric(rfstime)
  ) %>%
  na.omit()

## 建模不使用 pid
df_model <- df %>% select(-pid,-X)

## =========================================================
## 2) Train/Test split (固定一份 test 來做 IBS + calibration + KM)
## =========================================================
n <- nrow(df_model)
idx_test <- sample.int(n, size = floor(0.30 * n))
test_df  <- df_model[idx_test, ]
train_df <- df_model[-idx_test, ]

## =========================================================
## 3) Tools: OOB standardized CRPS + trapezoid IBS
##    CRPS/IBS 定義：CRPS = IBS / time (越小越好)
## =========================================================
get_oob_crps_std <- function(fit) {
  bs <- randomForestSRC:::get.brier.survival(fit)
  ## 不同版本欄位可能略不同；若報錯請先 str(bs) 看欄位
  if (!("crps.std" %in% names(bs))) stop("Cannot find crps.std in get.brier.survival output.")
  unname(bs$crps.std)
}

trapz <- function(x, y) {
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

## =========================================================
## 4) RSF coarse tuning (照你截圖範圍)
##    mtry: 1-8
##    nodesize: {5,10,15,20,30,50}
##    ntree: 1000
##    nsplit: {5,10}
##    指標：主 OOB std-CRPS；次 OOB C-index(=1-OOB error)
## =========================================================
grid1 <- expand.grid(
  mtry     = 1:8,
  nodesize = c(5, 10, 15, 20, 30, 50),
  nsplit   = c(5, 10),
  ntree    = 1000,
  stringsAsFactors = FALSE
)

fit_one_oob <- function(g) {
  fit <- rfsrc(
    Surv(rfstime, status) ~ .,
    data      = train_df,
    ntree     = g$ntree,
    mtry      = g$mtry,
    nodesize  = g$nodesize,
    nsplit    = g$nsplit,
    splitrule = "logrank",
    forest    = TRUE,
    importance = "none"
  )
  
  ## survival 的 OOB error 常用 PE = 1 - C
  oob_err    <- tail(fit$err.rate, 1)
  oob_cindex <- 1 - oob_err
  
  oob_crps_std <- get_oob_crps_std(fit)
  
  cbind(g, oob_err = oob_err, oob_cindex = oob_cindex, oob_crps_std = oob_crps_std)
}

cat("== Coarse tuning start ==\n")
res_list <- vector("list", nrow(grid1))
for (i in seq_len(nrow(grid1))) {
  g <- grid1[i, ]
  cat(sprintf("Grid %d/%d: mtry=%d nodesize=%d nsplit=%d\n", i, nrow(grid1), g$mtry, g$nodesize, g$nsplit))
  res_list[[i]] <- fit_one_oob(g)
}
res1 <- as.data.frame(do.call(rbind, res_list))

## rank-sum：CRPS 越小越好；C-index 越大越好
res1$rank_crps <- rank(res1$oob_crps_std, ties.method = "min")
res1$rank_c    <- rank(-res1$oob_cindex, ties.method = "min")
res1$rank_sum  <- res1$rank_crps + res1$rank_c

res1 <- res1[order(res1$rank_sum, res1$oob_crps_std, -res1$oob_cindex), ]
best_hp <- res1[1, c("mtry","nodesize","nsplit","ntree","oob_cindex","oob_crps_std")]
cat("\n== Best hyperparams (coarse) ==\n")
print(best_hp)

## =========================================================
## 5) Generate "fake measurement error" with 0.8
##    (A) Continuous: reliability r=0.8  ->  W = X + U, Var(U)=Var(X)*(1-r)/r
##    (B) Categorical (optional): grade misclassification tuned to Cohen’s kappa ≈ 0.8
## =========================================================

## (A) continuous ME generator
add_me_positive <- function(dat, var_names, r = 0.8,
                            method = c("log1p", "truncate"),
                            integer_vars = c("nodes")) {
  method <- match.arg(method)
  out <- dat
  
  for (vn in var_names) {
    x <- out[[vn]]
    v <- var(x, na.rm = TRUE)
    if (is.na(v) || v == 0) next
    
    if (method == "truncate") {
      # 原本做法 + 截斷到 0
      sd_u <- sqrt(v * (1 - r) / r)
      w <- x + rnorm(length(x), 0, sd_u)
      w <- pmax(0, w)
      
    } else {
      # 推薦：log1p 尺度加誤差，回到原尺度永遠 >= 0
      z <- log1p(x)
      vz <- var(z, na.rm = TRUE)
      sd_e <- sqrt(vz * (1 - r) / r)
      z_star <- z + rnorm(length(z), 0, sd_e)
      w <- exp(z_star) - 1
      w <- pmax(0, w)
    }
    
    # 若是計數型（nodes），做整數化
    if (vn %in% integer_vars) w <- as.integer(round(w))
    
    out[[vn]] <- w
  }
  
  out
}

## (B) categorical ME generator to target Cohen's kappa ≈ target_kappa
##     這裡用「只錯到相鄰級」的方式（對 grade 比較合理）
misclassify_grade_to_kappa <- function(x, target_kappa = 0.8, max_iter = 40, tol = 0.01) {
  x <- factor(x)
  lev <- levels(x)
  if (length(lev) < 2) return(list(x_star = x, kappa = 1, p = 0))
  
  ## 用二分搜尋找一個誤分類機率 p 使得 kappa 接近 target
  low <- 0; high <- 0.40  # 高一致性：通常不需要超過 0.4 的錯誤率
  best <- list(x_star = x, kappa = 1, p = 0, diff = Inf)
  
  for (iter in seq_len(max_iter)) {
    p <- (low + high) / 2
    x_star <- as.character(x)
    
    ## 隨機抽一部分做誤分類（只往相鄰級移動）
    flip <- runif(length(x_star)) < p
    for (i in which(flip)) {
      cur <- match(x_star[i], lev)
      if (cur == 1) {
        x_star[i] <- lev[2]
      } else if (cur == length(lev)) {
        x_star[i] <- lev[length(lev) - 1]
      } else {
        x_star[i] <- sample(c(lev[cur - 1], lev[cur + 1]), size = 1)
      }
    }
    
    x_star <- factor(x_star, levels = lev)
    
    k <- irr::kappa2(cbind(x, x_star))$value
    d <- abs(k - target_kappa)
    
    if (d < best$diff) best <- list(x_star = x_star, kappa = k, p = p, diff = d)
    if (d <= tol) break
    
    ## 如果目前 kappa 太高，代表誤差不夠 -> 增加 p；反之減少 p
    if (k > target_kappa) low <- p else high <- p
  }
  
  best[c("x_star","kappa","p")]
}

## 你要加誤差的變數（可按你研究需求調整）
cont_vars <- c("size", "pgr", "er", "nodes")  # 常見較會有量測誤差的連續欄位
use_grade_kappa <- TRUE                      # 要不要把 grade 做「kappa≈0.8」誤分類

## 產生一份「含測量誤差」的資料（整份資料都加；你也可只對 train 加）
set.seed(2026)
df_me <- add_me_positive(df_model, cont_vars, r = 0.8, method = "log1p")

if (use_grade_kappa) {
  set.seed(2026)
  k_out <- misclassify_grade_to_kappa(df_me$grade, target_kappa = 0.8, tol = 0.01)
  df_me$grade <- k_out$x_star
  cat(sprintf("\nGrade misclassification tuned: kappa=%.3f using p=%.3f\n", k_out$kappa, k_out$p))
}

## 用同一個切分索引，保持比較公平
train_me <- df_me[-idx_test, ]
test_me  <- df_me[idx_test, ]

## =========================================================
## 6) Fit RSF on ME data using best hyperparams
## =========================================================
best_fit <- rfsrc(
  Surv(rfstime, status) ~ .,
  data      = train_me,
  ntree     = best_hp$ntree,
  mtry      = best_hp$mtry,
  nodesize  = best_hp$nodesize,
  nsplit    = best_hp$nsplit,
  splitrule = "logrank",
  forest    = TRUE,
  importance = "permute"
)

print(best_fit)

## =========================================================
## 7) Evaluate: (Main) IBS/CRPS, (Second) C-index, (Must) calibration, + KM plot
## =========================================================
## 評估時間（例如 5 年內；依資料最大追蹤調整）
tau <- min(5 * 365, max(test_me$rfstime))
times_eval <- unique(pmax(30, round(seq(30, tau, length.out = 50))))

## riskRegression 會用 IPCW 方式處理 censoring，適合 IBS/Brier 類評估
sc <- Score(
  object  = list(RSF = best_fit),
  formula = Surv(rfstime, status) ~ 1,
  data    = test_me,
  times   = times_eval,
  metrics = c("brier", "auc"),
  plots   = "cal",
  cens.method = "ipcw"
)

## 取 Brier(t) 算 IBS（用梯形積分近似）
brier_df <- as.data.frame(sc$Brier$score)
brier_rsf <- subset(brier_df, model == "RSF")
ibs <- trapz(brier_rsf$times, brier_rsf$Brier) / (max(brier_rsf$times) - min(brier_rsf$times))
cat(sprintf("\n[Test] IBS (approx) = %.4f (lower is better)\n", ibs))

## C-index：用 RSF 預測 risk score（越大通常風險越高）
pred <- predict(best_fit, newdata = test_me)
risk_score <- pred$predicted
cobj <- survival::concordance(Surv(test_me$rfstime, test_me$status) ~ risk_score)
cat(sprintf("[Test] C-index = %.4f (higher is better)\n", cobj$concordance))

## Calibration curve（必做）
plotCalibration(sc, models = "RSF", times = max(times_eval), method = "quantile", q = 10)

## KM 圖：把人分成低/中/高風險組（例如三分位）
grp <- cut(risk_score,
           breaks = quantile(risk_score, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
           include.lowest = TRUE,
           labels = c("Low risk", "Mid risk", "High risk"))

km_fit <- survfit(Surv(rfstime, status) ~ grp, data = test_me)

plot(km_fit, col = 1:3, lwd = 2, xlab = "Days", ylab = "Survival probability",
     main = "Kaplan–Meier by RSF-predicted risk groups (ME data)")
legend("bottomleft", legend = levels(grp), col = 1:3, lwd = 2, bty = "n")

## =========================================================
## 8) (Optional next) Cox + SIMEX (measurement error correction)
##    你之後要做 Cox 校正比較時，把這段接上即可
## =========================================================
## install.packages("simex")
## library(simex)
## cox_naive <- coxph(Surv(rfstime, status) ~ age + meno + size + grade + nodes + pgr + er + hormon, data = train_me)
## simex_fit <- simex(cox_naive,
##                    SIMEXvariable = c("size","pgr","er","nodes"),
##                    measurement.error = sqrt(sapply(train_me[,c("size","pgr","er","nodes")], var) * (1-0.8)/0.8),
##                    lambda = seq(0, 2, by = 0.5),
##                    B = 100)
## summary(simex_fit)


