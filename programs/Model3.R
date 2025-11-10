##### Logistics Regression Model #####

load("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/data.RData")
data <- data_model

library(dplyr)
library(themis)   #step_smote()
library(tidymodels)
library(purrr)
library(yardstick)
library(glmnet)
library(recipes)

# =====================================================
# 1. Wrapper function for one run
# =====================================================
run_logistic_model <- function(seed, data, n_top = 20) {
  set.seed(seed)
  
  # -------------------
  # Train/test split
  # -------------------
  spl   <- initial_split(data, prop = 0.8, strata = ar_primary)
  train <- training(spl)
  test  <- testing(spl)
  
  # -------------------
  # Preprocessing (dummy first, then safe interactions)
  # -------------------
  rec <- recipe(ar_primary ~ ., data = train) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%   # dummy first
    step_interact(
      ~ starts_with("allergy_any_"):starts_with("race_") +
        starts_with("animal_allergy_"):starts_with("any_pet_") +
        starts_with("allergy_any_"):cot_ng_ml +
        LBXIGE:starts_with("race_")
    ) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_smote(ar_primary, over_ratio = 1.5, neighbors = 5)
  
  rec_prep <- prep(rec, training = train)
  train_baked <- bake(rec_prep, new_data = train)
  test_baked  <- bake(rec_prep, new_data = test)
  
  X_train <- as.matrix(train_baked %>% select(-ar_primary))
  y_train <- ifelse(train_baked$ar_primary == "1", 1, 0)
  X_test  <- as.matrix(test_baked %>% select(-ar_primary))
  y_test  <- ifelse(test_baked$ar_primary == "1", 1, 0)
  
  # -------------------
  # Penalty factor
  # -------------------
  pf <- rep(1, ncol(X_train))
  names(pf) <- colnames(X_train)
  pf[grep("allergy_any", names(pf))] <- 0
  pf[grep("race", names(pf))]        <- 0
  
  # -------------------
  # Fit penalized logistic regression
  # -------------------
  cv_fit <- cv.glmnet(
    X_train, y_train,
    family = "binomial",
    alpha = 0.5,
    penalty.factor = pf,
    nfolds = 5
  )
  
  # -------------------
  # Predictions + threshold tuning
  # -------------------
  pred_prob <- predict(cv_fit, newx = X_test, s = "lambda.min", type = "response")
  y_test_fac <- factor(y_test, levels = c(0,1))
  
  scan_threshold <- function(probs, truth, thresholds = seq(0.01, 0.99, 0.01)) {
    truth <- factor(truth, levels = c(0,1))
    map_dfr(thresholds, function(t) {
      pred_cls <- ifelse(probs >= t, 1, 0)
      eval_tbl <- tibble(
        truth = truth,
        pred_cls = factor(pred_cls, levels = c(0,1))
      )
      if (sum(eval_tbl$pred_cls == "1") == 0) {
        tibble(threshold = t, recall = NA, precision = NA, f1 = NA, acc = NA)
      } else {
        tibble(
          threshold = t,
          recall = sensitivity(eval_tbl, truth, pred_cls, event_level="second")$.estimate,
          precision = precision(eval_tbl, truth, pred_cls, event_level="second")$.estimate,
          f1 = f_meas(eval_tbl, truth, pred_cls, event_level="second")$.estimate,
          acc = accuracy(eval_tbl, truth, pred_cls)$.estimate
        )
      }
    })
  }
  
  met <- scan_threshold(pred_prob, y_test_fac)
  best_thr <- met %>%
    filter(recall >= 0.7) %>%
    slice_max(f1) %>%
    pull(threshold)
  
  pred_cls_thr <- ifelse(pred_prob >= best_thr, 1, 0)
  eval_tbl_thr <- tibble(
    truth = factor(y_test, levels = c(0,1)),
    .pred_1 = as.numeric(pred_prob),
    .pred_class = factor(pred_cls_thr, levels = c(0,1))
  )
  
  metrics <- tibble(
    seed = seed,
    acc  = accuracy(eval_tbl_thr, truth, .pred_class)$.estimate,
    ba   = bal_accuracy(eval_tbl_thr, truth, .pred_class)$.estimate,
    recall = sensitivity(eval_tbl_thr, truth, .pred_class, event_level="second")$.estimate,
    spec   = specificity(eval_tbl_thr, truth, .pred_class, event_level="second")$.estimate,
    precision = precision(eval_tbl_thr, truth, .pred_class, event_level="second")$.estimate,
    f1   = f_meas(eval_tbl_thr, truth, .pred_class, event_level="second")$.estimate,
    rocA = roc_auc(eval_tbl_thr, truth, .pred_1, event_level="second")$.estimate
  )
  
  # -------------------
  # Top N features
  # -------------------
  coef_table <- coef(cv_fit, s = "lambda.min")
  coef_table <- data.frame(
    variable = rownames(coef_table),
    coefficient = as.numeric(coef_table)
  ) %>%
    filter(variable != "(Intercept)", coefficient != 0) %>%
    arrange(desc(abs(coefficient)))
  
  top_features <- head(coef_table, n_top) %>%
    mutate(seed = seed)
  
  return(list(metrics = metrics, features = top_features))
}




# =====================================================
# 2. Run across your 10 random seeds
# =====================================================
# Example: suppose you already generated 10 random numbers from MD5
random_seeds <- c(44542169,162339537,296891153,2080800158,
                  825781625,6315874,1525767726,1579953055,
                  1794694103,1214191569)

# Run model across seeds
all_results <- map(random_seeds, ~ run_logistic_model(.x, data))

# Extract metrics
results <- bind_rows(map(all_results, "metrics"))

# Summarize metrics (mean/median/sd)
summary_results <- results %>%
  summarise(
    across(-seed, list(mean = mean, median = median, sd = sd),
           .names = "{.col}_{.fn}")
  )

# Extract top features across all seeds
features_across <- bind_rows(map(all_results, "features"))

# Count frequency of feature appearance
feature_rank_summary <- features_across %>%
  count(variable, sort = TRUE)

print(summary_results)        # averaged metrics
print(feature_rank_summary)   # most frequent top features

write.csv(summary_results, "result_log.csv")
