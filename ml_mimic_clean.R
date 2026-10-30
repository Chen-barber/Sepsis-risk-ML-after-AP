rm(list = ls())
gc()

SEED <- 123L
OUTCOME <- "Sepsis3"
OUT_DIR <- "data_ml_clean"
N_BOOT <- 1000L
THRESH_FIXED <- 0.5
THRESH <- THRESH_FIXED
N_THREADS <- max(1L, parallel::detectCores(logical = TRUE) - 2L)
Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OPENBLAS_NUM_THREADS = as.character(N_THREADS),
  MKL_NUM_THREADS = as.character(N_THREADS)
)

MODEL_ORDER <- c("KNN", "LR", "SLNN", "XGBoost", "LightGBM", "CatBoost", "GAMBoost")
feats <- c(
  "BUN_24h", "Calcium_24h", "Chloride_24h", "MV",
  "Norepinephrine", "PLT_24h", "SBP_24h"
)

KNN_K <- 201L
LR_LAMBDA <- 0.28503
SLNN_SIZE <- 37L
SLNN_DECAY <- 1
XGB_DEPTH <- 2L
XGB_ETA <- 0.03
XGB_NROUNDS <- 30L
XGB_LAMBDA <- 2
LGB_LEAVES <- 5L
LGB_LR <- 0.02
LGB_NROUNDS <- 15L
LGB_MIN_DATA <- 50L
CAT_DEPTH <- 3L
CAT_ITERS <- 40L
CAT_LR <- 0.05
CAT_L2 <- 5
GAMB_MSTOP <- 205L
GAMB_NU <- 0.09

local({
  r <- R.version
  pl <- file.path(
    "renv", "library", "windows",
    paste0("R-", r$major, ".", strsplit(r$minor, ".", fixed = TRUE)[[1]][1]),
    r$platform
  )
  if (dir.exists(pl)) .libPaths(c(normalizePath(pl), .libPaths()))
})

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(caret)
  library(nnet)
  library(kknn)
  library(glmnet)
  library(xgboost)
  library(lightgbm)
  library(mboost)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
})

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

read_xy <- function(path) {
  d <- read_csv(path, show_col_types = FALSE)
  d <- d[, c(feats, OUTCOME), drop = FALSE]
  d[feats] <- lapply(d[feats], as.numeric)
  d[[OUTCOME]] <- as.integer(d[[OUTCOME]])
  stopifnot(!anyNA(d))
  d
}

train <- read_xy(file.path("data_sets", "mimic_train.csv"))
test <- read_xy(file.path("data_sets", "mimic_internal_val.csv"))
ext <- read_xy(file.path("data_sets", "wmu_external_val.csv"))

make_fac <- function(d) {
  out <- d
  out[[OUTCOME]] <- factor(out[[OUTCOME]], levels = c(0L, 1L), labels = c("neg", "pos"))
  out
}
tr_f <- make_fac(train)
te_f <- make_fac(test)
ex_f <- make_fac(ext)
fml <- as.formula(paste(OUTCOME, "~", paste(feats, collapse = " + ")))

prob_pos <- function(pred_df) {
  if (is.data.frame(pred_df) && "pos" %in% names(pred_df)) return(as.numeric(pred_df$pos))
  if (is.matrix(pred_df) && "pos" %in% colnames(pred_df)) return(as.numeric(pred_df[, "pos"]))
  as.numeric(pred_df)
}

auc1 <- function(y, p) {
  y <- as.integer(y)
  if (length(unique(y)) < 2L) return(NA_real_)
  as.numeric(auc(roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<")))
}

auprc1 <- function(y, p) {
  y <- as.integer(y)
  if (sum(y == 1L) < 1L || sum(y == 0L) < 1L) return(NA_real_)
  as.numeric(pr.curve(
    scores.class0 = p[y == 1L], scores.class1 = p[y == 0L], curve = FALSE
  )$auc.integral)
}

cls_at_thresh <- function(y, p, thr = THRESH) {
  y <- as.integer(y)
  yhat <- as.integer(p >= thr)
  tp <- sum(yhat == 1L & y == 1L)
  tn <- sum(yhat == 0L & y == 0L)
  fp <- sum(yhat == 1L & y == 0L)
  fn <- sum(yhat == 0L & y == 1L)
  c(
    Accuracy = (tp + tn) / length(y),
    NPV = if ((tn + fn) == 0) NA_real_ else tn / (tn + fn),
    PPV = if ((tp + fp) == 0) NA_real_ else tp / (tp + fp),
    PRECISION = if ((tp + fp) == 0) NA_real_ else tp / (tp + fp),
    RECALL = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
    SENSITIVITY = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
    SPECIFICITY = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  )
}

auc_ci <- function(y, p) {
  a0 <- auc1(y, p)
  r <- tryCatch(
    roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<"),
    error = function(e) NULL
  )
  if (is.null(r)) return(c(AUC = a0, lo = NA_real_, hi = NA_real_))
  ci <- tryCatch(as.numeric(ci.auc(r, method = "delong")), error = function(e) NULL)
  if (is.null(ci) || length(ci) < 3L) return(c(AUC = a0, lo = NA_real_, hi = NA_real_))
  c(AUC = ci[2], lo = ci[1], hi = ci[3])
}

boot_metric_ci <- function(y, p, thr = THRESH, n_boot = N_BOOT, seed = SEED) {
  y <- as.integer(y)
  p <- as.numeric(p)
  point <- cls_at_thresh(y, p, thr)
  auc_v <- auc_ci(y, p)
  prauc0 <- auprc1(y, p)
  set.seed(seed)
  n <- length(y)
  boot_mat <- replicate(n_boot, {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[i])) < 2L) return(c(point, PRAUC = prauc0))
    c(cls_at_thresh(y[i], p[i], thr), PRAUC = auprc1(y[i], p[i]))
  })
  qlo <- apply(boot_mat, 1L, quantile, probs = 0.025, na.rm = TRUE)
  qhi <- apply(boot_mat, 1L, quantile, probs = 0.975, na.rm = TRUE)
  list(
    Accuracy = c(est = unname(point["Accuracy"]), lo = unname(qlo["Accuracy"]), hi = unname(qhi["Accuracy"])),
    AUC = c(est = unname(auc_v["AUC"]), lo = unname(auc_v["lo"]), hi = unname(auc_v["hi"])),
    NPV = c(est = unname(point["NPV"]), lo = unname(qlo["NPV"]), hi = unname(qhi["NPV"])),
    PPV = c(est = unname(point["PPV"]), lo = unname(qlo["PPV"]), hi = unname(qhi["PPV"])),
    PRAUC = c(est = prauc0, lo = unname(qlo["PRAUC"]), hi = unname(qhi["PRAUC"])),
    PRECISION = c(est = unname(point["PRECISION"]), lo = unname(qlo["PRECISION"]), hi = unname(qhi["PRECISION"])),
    RECALL = c(est = unname(point["RECALL"]), lo = unname(qlo["RECALL"]), hi = unname(qhi["RECALL"])),
    SENSITIVITY = c(est = unname(point["SENSITIVITY"]), lo = unname(qlo["SENSITIVITY"]), hi = unname(qhi["SENSITIVITY"])),
    SPECIFICITY = c(est = unname(point["SPECIFICITY"]), lo = unname(qlo["SPECIFICITY"]), hi = unname(qhi["SPECIFICITY"]))
  )
}

iso_fit_apply <- function(p_train, y_train, p_new) {
  o <- order(p_train)
  ir <- isoreg(p_train[o], y_train[o])
  xf <- ir$x[ir$iKnots]
  yf <- ir$yf[ir$iKnots]
  pmin(pmax(approx(xf, yf, xout = p_new, rule = 2)$y, 1e-6), 1 - 1e-6)
}

set.seed(SEED)
fit_knn <- train(
  fml, data = tr_f, method = "kknn",
  trControl = trainControl(method = "none", classProbs = TRUE),
  tuneGrid = data.frame(kmax = KNN_K, distance = 2, kernel = "rectangular"),
  preProcess = c("center", "scale")
)

x_tr <- scale(as.matrix(train[feats]))
x_te <- scale(
  as.matrix(test[feats]),
  center = attr(x_tr, "scaled:center"),
  scale = attr(x_tr, "scaled:scale")
)
x_ex <- scale(
  as.matrix(ext[feats]),
  center = attr(x_tr, "scaled:center"),
  scale = attr(x_tr, "scaled:scale")
)
set.seed(SEED)
fit_lr <- glmnet(
  x_tr, train[[OUTCOME]], family = "binomial",
  alpha = 0, lambda = LR_LAMBDA
)

set.seed(SEED)
fit_nn <- train(
  fml, data = tr_f, method = "nnet",
  trControl = trainControl(method = "none", classProbs = TRUE),
  tuneGrid = data.frame(size = SLNN_SIZE, decay = SLNN_DECAY),
  preProcess = c("center", "scale"),
  trace = FALSE, MaxNWts = 2000, maxit = 500
)

dtrain_xgb <- xgb.DMatrix(as.matrix(train[feats]), label = train[[OUTCOME]])
set.seed(SEED)
fit_xgb <- xgb.train(
  params = list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = XGB_DEPTH,
    eta = XGB_ETA,
    gamma = 0,
    lambda = XGB_LAMBDA,
    subsample = 0.8,
    colsample_bytree = 0.8,
    nthread = N_THREADS
  ),
  data = dtrain_xgb,
  nrounds = XGB_NROUNDS,
  verbose = 0
)

set.seed(SEED)
fit_lgb <- lgb.train(
  params = list(
    objective = "binary",
    metric = "auc",
    boosting = "gbdt",
    num_leaves = LGB_LEAVES,
    learning_rate = LGB_LR,
    feature_fraction = 0.95,
    bagging_fraction = 0.8,
    bagging_freq = 1L,
    min_data_in_leaf = LGB_MIN_DATA,
    lambda_l2 = 5,
    feature_pre_filter = FALSE,
    num_threads = N_THREADS,
    verbosity = -1L,
    seed = SEED
  ),
  data = lgb.Dataset(as.matrix(train[feats]), label = train[[OUTCOME]]),
  nrounds = LGB_NROUNDS
)

use_reticulate_cat <- !requireNamespace("catboost", quietly = TRUE)
if (use_reticulate_cat) {
  library(reticulate)
  cb <- import("catboost")
  set.seed(SEED)
  fit_cat <- cb$CatBoostClassifier(
    depth = as.integer(CAT_DEPTH),
    iterations = as.integer(CAT_ITERS),
    learning_rate = CAT_LR,
    l2_leaf_reg = CAT_L2,
    thread_count = N_THREADS,
    random_seed = SEED,
    verbose = FALSE,
    loss_function = "Logloss",
    eval_metric = "AUC",
    allow_writing_files = FALSE
  )
  fit_cat$fit(as.matrix(train[feats]), train[[OUTCOME]])
} else {
  set.seed(SEED)
  fit_cat <- catboost::catboost.train(
    catboost::catboost.load_pool(as.matrix(train[feats]), label = train[[OUTCOME]]),
    params = list(
      loss_function = "Logloss",
      eval_metric = "AUC",
      depth = CAT_DEPTH,
      iterations = CAT_ITERS,
      learning_rate = CAT_LR,
      l2_leaf_reg = CAT_L2,
      thread_count = N_THREADS,
      random_seed = SEED,
      logging_level = "Silent",
      allow_writing_files = FALSE
    )
  )
}

pred_cat <- function(d) {
  if (!use_reticulate_cat) {
    as.numeric(catboost::catboost.predict(
      fit_cat,
      catboost::catboost.load_pool(as.matrix(d[feats])),
      prediction_type = "Probability"
    ))
  } else {
    as.numeric(fit_cat$predict_proba(as.matrix(d[feats]))[, 2])
  }
}

gb_form <- as.formula(
  paste(OUTCOME, "~", paste(sprintf("btree(%s)", feats), collapse = " + "))
)
set.seed(SEED)
fit_gb <- gamboost(
  gb_form,
  data = tr_f,
  family = Binomial(),
  control = boost_control(mstop = GAMB_MSTOP, nu = GAMB_NU, trace = FALSE)
)

pred_raw <- function(model_name, d_num, d_fac) {
  switch(
    model_name,
    KNN = prob_pos(predict(fit_knn, newdata = d_fac, type = "prob")),
    LR = as.numeric(predict(
      fit_lr,
      scale(
        as.matrix(d_num[feats]),
        center = attr(x_tr, "scaled:center"),
        scale = attr(x_tr, "scaled:scale")
      ),
      type = "response"
    )),
    SLNN = prob_pos(predict(fit_nn, newdata = d_fac, type = "prob")),
    XGBoost = as.numeric(predict(fit_xgb, xgb.DMatrix(as.matrix(d_num[feats])))),
    LightGBM = as.numeric(predict(fit_lgb, as.matrix(d_num[feats]))),
    CatBoost = pred_cat(d_num),
    GAMBoost = as.numeric(predict(fit_gb, newdata = d_fac, type = "response")),
    stop(model_name)
  )
}

cohorts <- list(
  train = list(num = train, fac = tr_f, y = train[[OUTCOME]]),
  test = list(num = test, fac = te_f, y = test[[OUTCOME]]),
  external = list(num = ext, fac = ex_f, y = ext[[OUTCOME]])
)

raw <- lapply(MODEL_ORDER, function(m) {
  setNames(lapply(names(cohorts), function(cn) {
    pred_raw(m, cohorts[[cn]]$num, cohorts[[cn]]$fac)
  }), names(cohorts))
})
names(raw) <- MODEL_ORDER

cal <- lapply(MODEL_ORDER, function(m) {
  p_tr <- raw[[m]]$train
  setNames(lapply(names(cohorts), function(cn) {
    iso_fit_apply(p_tr, cohorts$train$y, raw[[m]][[cn]])
  }), names(cohorts))
})
names(cal) <- MODEL_ORDER

pred_long <- do.call(rbind, lapply(MODEL_ORDER, function(m) {
  do.call(rbind, lapply(names(cohorts), function(cn) {
    data.frame(
      model = m,
      cohort = cn,
      y_true = cohorts[[cn]]$y,
      y_prob_raw = raw[[m]][[cn]],
      y_prob = cal[[m]][[cn]],
      stringsAsFactors = FALSE
    )
  }))
}))
write_csv(pred_long, file.path(OUT_DIR, "predictions.csv"))

best_youden_thr <- function(y, p) {
  r <- roc(as.integer(y), as.numeric(p), quiet = TRUE, levels = c(0, 1), direction = "<")
  co <- coords(
    r, x = "best", best.method = "youden",
    ret = "threshold", transpose = FALSE
  )
  thr <- as.numeric(co$threshold[1])
  if (is.finite(thr)) return(thr)
  THRESH_FIXED
}
THRESH_YOUDEN <- best_youden_thr(cohorts$train$y, cal$GAMBoost$train)
THRESH <- THRESH_YOUDEN

compute_metrics_at_thr <- function(thr, thr_label) {
  metrics_list <- list()
  for (m in MODEL_ORDER) {
    for (cn in names(cohorts)) {
      y <- cohorts[[cn]]$y
      p <- cal[[m]][[cn]]
      ci_all <- boot_metric_ci(y, p, thr = thr)
      yhat <- as.integer(p >= thr)
      tp <- sum(yhat == 1L & y == 1L)
      fp <- sum(yhat == 1L & y == 0L)
      fn <- sum(yhat == 0L & y == 1L)
      f1 <- if ((2 * tp + fp + fn) == 0) 0 else 2 * tp / (2 * tp + fp + fn)
      metrics_list[[length(metrics_list) + 1L]] <- data.frame(
        threshold_rule = thr_label,
        threshold_value = thr,
        model = m,
        cohort = cn,
        n = length(y),
        event_rate = mean(y),
        Accuracy = ci_all$Accuracy["est"],
        Accuracy_CI_low = ci_all$Accuracy["lo"],
        Accuracy_CI_high = ci_all$Accuracy["hi"],
        AUC = ci_all$AUC["est"],
        AUC_CI_low = ci_all$AUC["lo"],
        AUC_CI_high = ci_all$AUC["hi"],
        NPV = ci_all$NPV["est"],
        NPV_CI_low = ci_all$NPV["lo"],
        NPV_CI_high = ci_all$NPV["hi"],
        PPV = ci_all$PPV["est"],
        PPV_CI_low = ci_all$PPV["lo"],
        PPV_CI_high = ci_all$PPV["hi"],
        PRAUC = ci_all$PRAUC["est"],
        PRAUC_CI_low = ci_all$PRAUC["lo"],
        PRAUC_CI_high = ci_all$PRAUC["hi"],
        F1 = f1,
        Brier = mean((p - y)^2),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(metrics_list)
}
metrics_fixed <- compute_metrics_at_thr(THRESH_FIXED, "thresh0.5")
metrics_youden <- compute_metrics_at_thr(THRESH_YOUDEN, "Youden")
metrics <- metrics_youden

params <- data.frame(
  model = MODEL_ORDER,
  param = c(
    sprintf("k=%d,distance=2,kernel=rectangular,scale=TRUE", KNN_K),
    sprintf("ridge glmnet alpha=0 lambda=%s", LR_LAMBDA),
    sprintf("size=%d,decay=%s", SLNN_SIZE, SLNN_DECAY),
    sprintf(
      "gbtree,max_depth=%d,eta=%s,nrounds=%d,gamma=0,lambda=%s,subsample=0.8,colsample=0.8,nthread=%d",
      XGB_DEPTH, XGB_ETA, XGB_NROUNDS, XGB_LAMBDA, N_THREADS
    ),
    sprintf(
      "gbdt,num_leaves=%d,lr=%s,nrounds=%d,feature_fraction=0.95,bagging=0.8,min_data=%d,lambda_l2=5,nthread=%d",
      LGB_LEAVES, LGB_LR, LGB_NROUNDS, LGB_MIN_DATA, N_THREADS
    ),
    sprintf(
      "depth=%d,iterations=%d,lr=%s,l2=%s,thread_count=%d",
      CAT_DEPTH, CAT_ITERS, CAT_LR, CAT_L2, N_THREADS
    ),
    sprintf("mstop=%d,nu=%s,btree,Binomial", GAMB_MSTOP, GAMB_NU)
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "metrics_thresh0.5")
writeData(wb, "metrics_thresh0.5", metrics_fixed)
addWorksheet(wb, "metrics_Youden")
writeData(wb, "metrics_Youden", metrics_youden)
addWorksheet(wb, "thresholds")
writeData(
  wb, "thresholds",
  data.frame(
    rule = c("thresh0.5", "Youden"),
    value = c(THRESH_FIXED, THRESH_YOUDEN),
    stringsAsFactors = FALSE
  )
)
addWorksheet(wb, "best_params")
writeData(wb, "best_params", params)
addWorksheet(wb, "features")
writeData(wb, "features", data.frame(variable = feats))
saveWorkbook(wb, file.path(OUT_DIR, "metrics.xlsx"), overwrite = TRUE)

MODEL_COLS <- c(
  KNN = "#4DBBD5", LR = "#E64B35", SLNN = "#00A087",
  XGBoost = "#3C5488", LightGBM = "#F39B7F", CatBoost = "#8491B4",
  GAMBoost = "#91D1C2"
)

make_roc_panel <- function(cn) {
  y <- cohorts[[cn]]$y
  curves <- do.call(rbind, lapply(MODEL_ORDER, function(m) {
    r <- roc(y, cal[[m]][[cn]], quiet = TRUE, levels = c(0, 1), direction = "<")
    data.frame(
      model = m,
      fpr = 1 - r$specificities,
      tpr = r$sensitivities,
      stringsAsFactors = FALSE
    )
  }))
  curves$model <- factor(curves$model, levels = MODEL_ORDER)
  ggplot(curves, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.55) +
    scale_colour_manual(values = MODEL_COLS) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = cn) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), legend.title = element_blank())
}

p_roc <- make_roc_panel("train") + make_roc_panel("test") + make_roc_panel("external")
ggsave(file.path(OUT_DIR, "roc_panels.pdf"), p_roc, width = 10.5, height = 3.4, device = pdf)

saveRDS(
  list(
    fit_knn = fit_knn,
    fit_lr = fit_lr,
    fit_nn = fit_nn,
    fit_xgb = fit_xgb,
    fit_lgb = fit_lgb,
    fit_cat = if (!use_reticulate_cat) fit_cat else NULL,
    fit_gb = fit_gb,
    feats = feats,
    params = params,
    x_tr_center = attr(x_tr, "scaled:center"),
    x_tr_scale = attr(x_tr, "scaled:scale")
  ),
  file.path(OUT_DIR, "models_bundle.rds")
)
