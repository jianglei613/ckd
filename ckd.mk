```{r}
required_pkgs <- c("tidyverse", "readxl", "caret", "ranger", "xgboost", 
                   "pROC", "ggplot2", "ROSE", "recipes")
for(pkg in required_pkgs) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
set.seed(123)

# 1. 读取数据 -------------------------------------------------
read_medical_data <- function(file_path, year) {
  data <- readxl::read_excel(file_path) %>%
    mutate(year = year) %>%
    rename_with(~ tolower(gsub(" ", "_", .x)))
  return(data)
}

data_2015 <- read_medical_data("15_data.xlsx", 2015)
data_2018 <- read_medical_data("18_data.xlsx", 2018) 
data_2020 <- read_medical_data("20_data.xlsx", 2020)

# 2. 初步清洗 ------------------------------------------------
continuous_vars <- c("albumin", "glu", "mcv", "ua", "globulin", "mch", "tc", 
                     "rdw_sd", "ldl", "hemoglobin", "total_protein", "neut_percent",
                     "rbc", "plt", "hdl_c", "tg", "bun", "creatinine", "mchc", 
                     "hct", "egfr")

clean_medical_data <- function(data) {
  data %>%
    distinct(id_card, year, .keep_all = TRUE) %>%
    mutate(male = ifelse(gender == 1, 1, ifelse(gender == 2, 0, NA))) %>%
    group_by(id_card) %>%
    filter(sum(is.na(cur_data())) / (ncol(.) * n()) < 0.5) %>%
    ungroup()
}

data_2015_clean <- clean_medical_data(data_2015)
data_2018_clean <- clean_medical_data(data_2018)
data_2020_clean <- clean_medical_data(data_2020)

# 3. 定义基线CKD并排除 -----------------------------------------
is_ckd <- function(egfr, urine_protein) {
  as.integer(egfr < 60 | urine_protein > 1)
}

baseline_ckd_ids <- data_2015_clean %>%
  mutate(ckd_baseline = is_ckd(egfr, urine_protein)) %>%
  filter(ckd_baseline == 1) %>%
  pull(id_card)

eligible_ids <- setdiff(unique(data_2015_clean$id_card), baseline_ckd_ids)

data_2015_eligible <- data_2015_clean %>% filter(id_card %in% eligible_ids)
data_2018_eligible <- data_2018_clean %>% filter(id_card %in% eligible_ids)
data_2020_eligible <- data_2020_clean %>% filter(id_card %in% eligible_ids)

# 4. 结局定义 ------------------------------------------------
incident_ckd_df <- data_2020_eligible %>%
  mutate(incident_ckd = is_ckd(egfr, urine_protein)) %>%
  select(id_card, incident_ckd) %>%
  mutate(incident_ckd = factor(incident_ckd, levels = c(0,1), labels = c("No","Yes")))

# 5. 特征工程（只使用2015和2018）-----------------------------
baseline_2015 <- data_2015_eligible %>%
  select(id_card, age, male, all_of(continuous_vars)) %>%
  rename_with(~ paste0(., "_2015"), .cols = all_of(continuous_vars))

data_2018_selected <- data_2018_eligible %>%
  select(id_card, all_of(continuous_vars)) %>%
  rename_with(~ paste0(., "_2018"), .cols = all_of(continuous_vars))

wide_data <- baseline_2015 %>%
  inner_join(data_2018_selected, by = "id_card") %>%
  inner_join(incident_ckd_df, by = "id_card")

feature_data <- wide_data %>%
  mutate(
    egfr_change_abs = egfr_2018 - egfr_2015,
    creatinine_change_abs = creatinine_2018 - creatinine_2015,
    egfr_change_rel = ifelse(egfr_2015 != 0, (egfr_2018 - egfr_2015) / egfr_2015, 0),
    creatinine_change_rel = ifelse(creatinine_2015 != 0, (creatinine_2018 - creatinine_2015) / creatinine_2015, 0),
    egfr_change_annual = egfr_change_abs / 3,
    creatinine_change_annual = creatinine_change_abs / 3,
    albumin_change = ifelse(albumin_2015 != 0, (albumin_2018 - albumin_2015) / albumin_2015, 0),
    bun_change = ifelse(bun_2015 != 0, (bun_2018 - bun_2015) / bun_2015, 0),
    albumin_creatinine_ratio_2015 = albumin_2015 / creatinine_2015,
    bun_creatinine_ratio_2015 = bun_2015 / creatinine_2015,
    baseline_egfr_category = case_when(
      egfr_2015 >= 90 ~ 0,
      egfr_2015 >= 60 ~ 1,
      TRUE ~ 2
    ),
    egfr_trend = case_when(
      egfr_change_annual < -3 ~ 2,
      egfr_change_annual < 0 ~ 1,
      TRUE ~ 0
    )
  )

model_data <- feature_data %>%
  select(-id_card) %>%
  rename(ckd_status = incident_ckd)

model_data$ckd_status <- factor(model_data$ckd_status, levels = c("No","Yes"))
cat("总样本量:", nrow(model_data), "，新发CKD率:", round(mean(model_data$ckd_status=="Yes")*100,2), "%\n")

# 6. 划分训练、验证、测试集（分层）-----------------------------
set.seed(123)
train_idx <- caret::createDataPartition(model_data$ckd_status, p = 0.7, list = FALSE)
temp <- model_data[train_idx, ]
test  <- model_data[-train_idx, ]

set.seed(123)
val_idx <- caret::createDataPartition(temp$ckd_status, p = 0.15/0.7, list = FALSE)
train <- temp[-val_idx, ]
val   <- temp[val_idx, ]

cat("训练集:", nrow(train), "正例:", sum(train$ckd_status=="Yes"),
    "\n验证集:", nrow(val), "正例:", sum(val$ckd_status=="Yes"),
    "\n测试集:", nrow(test), "正例:", sum(test$ckd_status=="Yes"), "\n")

# 7. 预处理（基于训练集，应用到验证/测试）-----------------------
predictors <- setdiff(names(train), "ckd_status")
num_pred <- predictors[sapply(train[predictors], is.numeric)]

rec <- recipes::recipe(ckd_status ~ ., data = train) %>%
  recipes::step_zv(all_predictors()) %>%
  recipes::step_impute_median(all_numeric_predictors()) %>%
  recipes::step_normalize(all_numeric_predictors())

rec_prep <- recipes::prep(rec, training = train)
train_processed <- recipes::bake(rec_prep, new_data = train)
val_processed   <- recipes::bake(rec_prep, new_data = val)
test_processed  <- recipes::bake(rec_prep, new_data = test)

# 移除零方差变量
remove_zero_var <- function(df) {
  vars <- names(df)[sapply(df, function(x) is.numeric(x) && length(unique(x)) == 1)]
  if(length(vars) > 0) {
    cat("移除零方差变量:", paste(vars, collapse=", "), "\n")
    df <- df %>% select(-all_of(vars))
  }
  return(df)
}
train_processed <- remove_zero_var(train_processed)
val_processed <- val_processed %>% select(all_of(names(train_processed)))
test_processed <- test_processed %>% select(all_of(names(train_processed)))

# 8. 处理类别不平衡（仅训练集，使用ROSE）-----------------------
if(mean(train_processed$ckd_status=="Yes") < 0.3) {
  set.seed(123)
  train_balanced <- ROSE::ROSE(ckd_status ~ ., data = train_processed, seed = 123)$data
  train_balanced$ckd_status <- factor(train_balanced$ckd_status, levels = c("No","Yes"))
  cat("ROSE后训练集样本量:", nrow(train_balanced), "正例:", sum(train_balanced$ckd_status=="Yes"), "\n")
} else {
  train_balanced <- train_processed
}

```
```{r}
# 9. 模型训练（手动调参，基于验证集AUC）------------------------
eval_val_auc <- function(pred_prob, actual) {
  roc_obj <- pROC::roc(actual, pred_prob, quiet = TRUE)
  return(as.numeric(pROC::auc(roc_obj)))
}

# 9.1 随机森林（ranger）
rf_param_grid <- c(2,4,6,8)
best_rf_auc <- 0
best_rf_mtry <- 2
best_rf_model <- NULL

for(mtry_val in rf_param_grid) {
  set.seed(123)
  rf_tmp <- ranger::ranger(
    ckd_status ~ ., data = train_balanced,
    num.trees = 500,
    mtry = mtry_val,
    importance = "impurity",
    probability = TRUE,
    seed = 123
  )
  prob_val <- predict(rf_tmp, val_processed, type = "response")$predictions[, "Yes"]
  auc_val <- eval_val_auc(prob_val, val_processed$ckd_status)
  cat("RF mtry =", mtry_val, "验证集AUC =", round(auc_val,4), "\n")
  if(auc_val > best_rf_auc) {
    best_rf_auc <- auc_val
    best_rf_mtry <- mtry_val
    best_rf_model <- rf_tmp
  }
}
cat("最佳RF: mtry =", best_rf_mtry, "验证集AUC =", round(best_rf_auc,4), "\n")

# 9.2 XGBoost
xgb_param_grid <- expand.grid(
  max_depth = c(3,5,7),
  eta = c(0.01, 0.05, 0.1)
)
best_xgb_auc <- 0
best_xgb_params <- NULL
best_xgb_model <- NULL

dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_balanced %>% select(-ckd_status)),
                               label = as.numeric(train_balanced$ckd_status == "Yes"))
dval   <- xgboost::xgb.DMatrix(data = as.matrix(val_processed %>% select(-ckd_status)),
                               label = as.numeric(val_processed$ckd_status == "Yes"))

for(i in 1:nrow(xgb_param_grid)) {
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = xgb_param_grid$max_depth[i],
    eta = xgb_param_grid$eta[i],
    verbose = 0
  )
  set.seed(123)
  xgb_tmp <- xgboost::xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
  pred_val <- predict(xgb_tmp, dval)
  auc_val <- eval_val_auc(pred_val, val_processed$ckd_status)
  cat("XGB depth=", params$max_depth, "eta=", params$eta, "验证集AUC=", round(auc_val,4), "\n")
  if(auc_val > best_xgb_auc) {
    best_xgb_auc <- auc_val
    best_xgb_params <- params
    best_xgb_model <- xgb_tmp
  }
}
cat("最佳XGB: depth=", best_xgb_params$max_depth, "eta=", best_xgb_params$eta, "AUC=", round(best_xgb_auc,4), "\n")

# 9.3 逻辑回归
lr_model <- glm(ckd_status ~ ., data = train_balanced, family = binomial())
prob_val_lr <- predict(lr_model, val_processed, type = "response")
auc_lr_val <- eval_val_auc(prob_val_lr, val_processed$ckd_status)
cat("逻辑回归验证集AUC =", round(auc_lr_val,4), "\n")

# 10. 在验证集上选择最优阈值（Youden指数）---------------------
get_best_threshold <- function(prob, actual) {
  roc_obj <- pROC::roc(actual, prob, quiet = TRUE)
  coords(roc_obj, "best", ret = "threshold")$threshold
}

thresh_rf <- get_best_threshold(
  predict(best_rf_model, val_processed, type = "response")$predictions[, "Yes"],
  val_processed$ckd_status
)
thresh_xgb <- get_best_threshold(
  predict(best_xgb_model, dval),
  val_processed$ckd_status
)
thresh_lr <- get_best_threshold(prob_val_lr, val_processed$ckd_status)
```
```{r}
# 11. 测试集最终评估 -------------------------------------------
evaluate_test <- function(prob_test, actual, threshold, model_name) {
  pred_class <- factor(ifelse(prob_test > threshold, "Yes", "No"), levels = c("No","Yes"))
  cm <- caret::confusionMatrix(pred_class, actual, positive = "Yes")
  roc_obj <- pROC::roc(actual, prob_test, quiet = TRUE)
  data.frame(
    Model = model_name,
    AUC = round(as.numeric(pROC::auc(roc_obj)), 3),
    Sensitivity = round(cm$byClass["Sensitivity"], 3),
    Specificity = round(cm$byClass["Specificity"], 3),
    Accuracy = round(cm$overall["Accuracy"], 3),
    F1 = round(cm$byClass["F1"], 3),
    PPV = round(cm$byClass["Pos Pred Value"], 3),
    NPV = round(cm$byClass["Neg Pred Value"], 3),
    Threshold = round(threshold, 3)
  )
}

prob_test_rf <- predict(best_rf_model, test_processed, type = "response")$predictions[, "Yes"]
prob_test_xgb <- predict(best_xgb_model, xgboost::xgb.DMatrix(data = as.matrix(test_processed %>% select(-ckd_status))))
prob_test_lr <- predict(lr_model, test_processed, type = "response")

results_rf <- evaluate_test(prob_test_rf, test_processed$ckd_status, thresh_rf, "Random Forest")
results_xgb <- evaluate_test(prob_test_xgb, test_processed$ckd_status, thresh_xgb, "XGBoost")
results_lr <- evaluate_test(prob_test_lr, test_processed$ckd_status, thresh_lr, "Logistic Regression")

final_perf <- dplyr::bind_rows(results_rf, results_xgb, results_lr)
print(final_perf)
```
```{r}
# 12. 可视化 ---------------------------------------------------
roc_rf <- pROC::roc(test_processed$ckd_status, prob_test_rf, quiet = TRUE)
roc_xgb <- pROC::roc(test_processed$ckd_status, prob_test_xgb, quiet = TRUE)
roc_lr <- pROC::roc(test_processed$ckd_status, prob_test_lr, quiet = TRUE)

# 修正：使用 pROC::ggroc，不是 ggplot2::ggroc
roc_plot <- pROC::ggroc(list(RF = roc_rf, XGB = roc_xgb, LR = roc_lr), size = 1) +
  ggplot2::geom_abline(linetype = "dashed", color = "gray") +
  ggplot2::theme_minimal() +
  #ggplot2::labs(title = "ROC Curve Comparison") +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave("ROC_Curve_Final.png", roc_plot, width = 8, height = 6, dpi = 300)

# 随机森林特征重要性
rf_imp <- ranger::importance(best_rf_model)
imp_df <- data.frame(Feature = names(rf_imp), Importance = rf_imp) %>%
  arrange(desc(Importance)) %>%
  head(15)
imp_plot <- ggplot2::ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() +
  #ggplot2::labs(x = "Feature", y = "Importance", title = "Feature Importance Ranking") +
  ggplot2::labs(x = "Feature", y = "Importance") +
  ggplot2::theme_minimal()
ggplot2::ggsave("Feature_Importance_Final.png", imp_plot, width = 10, height = 7, dpi = 300)


# 13. 保存结果 -----------------------------------------------
saveRDS(best_rf_model, "RF_model_final.rds")
saveRDS(best_xgb_model, "XGB_model_final.rds")
saveRDS(lr_model, "LR_model_final.rds")
write.csv(final_perf, "Performance_metrics_final.csv", row.names = FALSE)

cat("\n=== 分析完成 ===\n")
cat("最终测试集性能:\n")
print(final_perf)
```

```{r}
# 风险分层（基于RF）
test_with_risk <- test_processed %>%
  mutate(
    pred_prob = prob_test_rf,
    risk_stratum = case_when(
      pred_prob < 0.1 ~ "Low (<10%)",
      pred_prob < 0.2 ~ "Medium (10-20%)",
      pred_prob < 0.3 ~ "High (20-30%)",
      TRUE ~ "Very high (>30%)"
    )
  )

# 将 risk_stratum 转换为有序因子（指定顺序）
test_with_risk <- test_with_risk %>%
  mutate(risk_stratum = factor(risk_stratum,
                               levels = c("Low (<10%)", "Medium (10-20%)",
                                          "High (20-30%)", "Very high (>30%)")))
risk_summary <- test_with_risk %>%
  group_by(risk_stratum) %>%
  summarise(n = n(), cases = sum(ckd_status == "Yes"), incidence = cases/n*100)
risk_plot <- ggplot2::ggplot(risk_summary, aes(x = risk_stratum, y = incidence, fill = risk_stratum)) +
  ggplot2::geom_col() + 
  ggplot2::geom_text(aes(label = paste0("n=", n, "\n", round(incidence,1),"%")), vjust = 0.5) +
  ggplot2::labs( #title = "Incidence of new-onset CKD in different risk strata",
      x = "Risk stratification",
      y = "Incidence rate of new CKD (%)",
      fill = "Risk stratification") +
  ggplot2::theme_minimal()
risk_plot

#ggplot2::ggsave("Risk_Stratification_Final.png", risk_plot, width = 6, height = 8, dpi = 300)

```


```{r}
# ============================================================
# 修正：生成三模型的校准曲线（Calibration Curve）
# ============================================================

# 准备三个模型的预测概率（假设已存在）
# prob_test_rf, prob_test_xgb, prob_test_lr 已在主代码中定义
# test_processed$ckd_status 为因子，需要转换为0/1数值

actual <- as.numeric(test_processed$ckd_status == "Yes")

# 定义函数：计算单个模型的校准曲线数据
get_calibration_data <- function(pred_prob, actual, model_name) {
  data.frame(Pred = pred_prob, Obs = actual) %>%
    mutate(bin = cut(Pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(Pred),
      mean_obs = mean(Obs),
      n = n(),
      Model = model_name,
      .groups = "drop"
    )
}

# 分别计算三个模型的校准数据
cal_rf <- get_calibration_data(prob_test_rf, actual, "Random Forest")
cal_xgb <- get_calibration_data(prob_test_xgb, actual, "XGBoost")
cal_lr  <- get_calibration_data(prob_test_lr, actual, "Logistic Regression")

cal_all <- bind_rows(cal_rf, cal_xgb, cal_lr)

# 绘制三条校准曲线
cal_plot <- ggplot(cal_all, aes(x = mean_pred, y = mean_obs, color = Model)) +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray30", size = 0.8) +
  scale_color_manual(values = c("Random Forest" = "#E41A1C", 
                                "XGBoost" = "#377EB8", 
                                "Logistic Regression" = "#4DAF4A")) +
  labs(x = "Predicted Probability", 
       y = "Observed Proportion",
       title = "Calibration Curves") +
  coord_equal() + 
  xlim(0, 0.5) + 
  ylim(0, 0.5) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("Calibration_Curve_Three_Models.png", cal_plot, width = 7, height = 6, dpi = 300)

# ============================================================
# 修正：生成三模型的决策曲线（包含 Treat All / Treat None）
# ============================================================

thresholds <- seq(0.01, 0.3, by = 0.01)

# 计算单个模型的净收益
calc_nb <- function(pred_prob, actual, thresholds) {
  sapply(thresholds, function(t) {
    mean((pred_prob > t) * actual) - mean((pred_prob > t) * (1 - actual)) * (t / (1 - t))
  })
}

# 计算参考线：Treat All（所有个体干预）和 Treat None（不干预）
treat_all_nb <- sapply(thresholds, function(t) mean(actual) - (1 - mean(actual)) * (t / (1 - t)))
treat_none_nb <- rep(0, length(thresholds))

# 计算三个模型的净收益
nb_rf  <- calc_nb(prob_test_rf, actual, thresholds)
nb_xgb <- calc_nb(prob_test_xgb, actual, thresholds)
nb_lr  <- calc_nb(prob_test_lr, actual, thresholds)

# 构建数据框
dc_data <- data.frame(
  Threshold = rep(thresholds, 5),
  NetBenefit = c(treat_all_nb, treat_none_nb, nb_rf, nb_xgb, nb_lr),
  Strategy = rep(c("Treat All", "Treat None", "Random Forest", "XGBoost", "Logistic Regression"), 
                 each = length(thresholds))
)

# 避免 Treat None 为 NA（保持为0）
dc_data$NetBenefit[is.na(dc_data$NetBenefit)] <- 0

# 绘制决策曲线
dc_plot <- ggplot(dc_data, aes(x = Threshold, y = NetBenefit, color = Strategy)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = c("Treat All" = "gray50", 
                                "Treat None" = "black",
                                "Random Forest" = "#E41A1C",
                                "XGBoost" = "#377EB8",
                                "Logistic Regression" = "#4DAF4A")) +
  labs(x = "Threshold Probability", 
       y = "Net Benefit") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ylim(-0.05, 0.2)

ggsave("Decision_Curve_Three_Models.png", dc_plot, width = 7, height = 6, dpi = 300)
```


```{r}
# 校准曲线
cal_data <- data.frame(
  Pred = prob_test_rf,
  Obs = as.numeric(test_processed$ckd_status == "Yes")
)
cal_data$bin <- cut(cal_data$Pred, breaks = seq(0,1,0.1), include.lowest = TRUE)
cal_sum <- cal_data %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(Pred), mean_obs = mean(Obs), n = n())
cal_plot <- ggplot(cal_sum, aes(x = mean_pred, y = mean_obs)) +
  geom_point(aes(size = n)) + geom_abline(linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Predicted Probability", y = "Observed Proportion") +
  coord_equal() + xlim(0,0.5) + ylim(0,0.5) +
  theme_minimal()
ggsave("Calibration_Curve.png", cal_plot, width = 6, height = 5)

# 决策曲线（简化版）
thresholds <- seq(0.01, 0.3, by = 0.01)
net_benefit <- sapply(thresholds, function(t) {
  mean((prob_test_rf > t) * (test_processed$ckd_status == "Yes")) - 
    mean((prob_test_rf > t) * (test_processed$ckd_status == "No")) * (t/(1-t))
})
dc_data <- data.frame(Threshold = thresholds, NetBenefit = net_benefit)
dc_plot <- ggplot(dc_data, aes(x = Threshold, y = NetBenefit)) +
  geom_line() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Net Benefit") + theme_minimal()
ggsave("Decision_Curve.png", dc_plot, width = 6, height = 5)

```


