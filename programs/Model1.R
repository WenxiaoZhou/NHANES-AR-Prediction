##################### Setup #####################
library(tidymodels)
library(tidyverse)
library(doFuture)
library(finetune)
library(themis)
library(ranger)
library(xgboost)
library(fastshap)
library(shapviz)
library(dplyr)
library(ggplot2)
library(htmltools)

# Parallel backend
n_cores <- parallel::detectCores(logical = TRUE)
registerDoFuture()
plan(multisession, workers = max(1, n_cores - 1))

# Load data
load("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/data.RData")
data <- data_model
data$ar_primary <- factor(data$ar_primary, levels = c(0, 1))


##################### Recipe #####################
rec <- recipe(ar_primary ~ ., data = data) %>%
  step_string2factor(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_other(all_nominal_predictors(), threshold = 0.03) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(~ starts_with("allergy_any"):starts_with("race") +
                  starts_with("animal_allergy"):starts_with("any_pet") +
                  starts_with("allergy_any"):cot_ng_ml +
                  LBXIGE:starts_with("race")) %>%
  step_log(LBXEOPCT, nlr, LBXIGE, crp_mg_l, cot_ng_ml, offset = 1) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_predictors()) %>%
  step_downsample(ar_primary, under_ratio = 1.5) %>%
  step_smote(ar_primary, over_ratio = 1.0, neighbors = 3)


##################### Model Specs #####################
# Random Forest
rf_tune <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_engine("ranger",
             probability = TRUE,
             importance = "permutation",
             num.threads = max(1, parallel::detectCores(logical = TRUE) - 1)) %>%
  set_mode("classification")

# XGBoost
boost_tune <- boost_tree(
  trees = 2000,
  tree_depth = tune(),
  learn_rate = tune(),
  min_n = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  mtry = tune()
) %>%
  set_engine("xgboost",
             nthread = max(1, parallel::detectCores(logical = TRUE) - 1)) %>%
  set_mode("classification")


##################### Helper Function #####################
best_thr_recall_f1 <- function(truth, scores, min_recall = 0.70,
                               ths = seq(0.01, 0.99, 0.01)) {
  cand <- purrr::map_dfr(ths, function(t) {
    cls <- factor(ifelse(scores >= t, "1", "0"), levels = c("0", "1"))
    tibble::tibble(
      threshold = t,
      recall = yardstick::sens_vec(truth, cls, event_level = "second"),
      ppv = yardstick::precision_vec(truth, cls, event_level = "second"),
      f1 = yardstick::f_meas_vec(truth, cls, event_level = "second", beta = 1)
    )
  })
  cand_main <- dplyr::filter(cand, recall >= min_recall)
  if (nrow(cand_main) > 0) return(cand_main$threshold[which.max(cand_main$f1)])
  best_recall <- max(cand$recall, na.rm = TRUE)
  cand_fallback <- dplyr::filter(cand, recall == best_recall)
  return(cand_fallback$threshold[which.max(cand_fallback$f1)])
}


##################### Evaluation Function #####################
evaluate_model_with_seed <- function(model_spec, model_name, seed, data, rec) {
  set.seed(seed)
  
  # Split
  spl   <- initial_split(data, prop = 0.8, strata = "ar_primary")
  train <- training(spl)
  test  <- testing(spl)
  
  # Workflow
  wf <- workflow() %>% add_recipe(rec) %>% add_model(model_spec)
  params <- wf %>% extract_parameter_set_dials()
  if ("mtry" %in% params$id) {
    params <- params %>%
      update(mtry = finalize(mtry(), train %>% select(-ar_primary)))
  }
  
  # Bayesian tuning
  folds <- vfold_cv(train, v = 5, strata = ar_primary)
  res <- tune_bayes(
    wf, resamples = folds, param_info = params,
    initial = 10, iter = 20,
    metrics = metric_set(bal_accuracy),
    control = control_bayes(no_improve = 10, verbose = FALSE, parallel_over = "resamples")
  )
  
  best <- select_best(res, metric = "bal_accuracy")
  final_wf <- finalize_workflow(wf, best)
  fit_final <- fit(final_wf, data = train)
  
  # Calibration
  preds_cal <- predict(fit_final, train, type = "prob") %>%
    bind_cols(train %>% select(ar_primary)) %>%
    transmute(ar_primary, score = .pred_1)
  cal_model <- logistic_reg() %>% set_engine("glm") %>% set_mode("classification")
  cal_fit   <- fit(cal_model, ar_primary ~ score, data = preds_cal)
  
  calibrate_preds <- function(fit_rf, fit_cal, newdata){
    rf_preds <- predict(fit_rf, newdata, type = "prob") %>%
      transmute(score = .pred_1)
    predict(fit_cal, rf_preds, type = "prob")$.pred_1
  }
  
  # Threshold
  val_res <- fit_resamples(final_wf, folds, control = control_resamples(save_pred = TRUE))
  preds_cv <- collect_predictions(val_res) %>%
    transmute(truth = fct_relevel(ar_primary, "0","1"), score = .pred_1)
  preds_cv$cal_score <- predict(cal_fit, new_data = preds_cv %>% select(score), type = "prob")$.pred_1
  thr <- best_thr_recall_f1(preds_cv$truth, preds_cv$cal_score, min_recall = 0.70)
  
  # Test set
  preds_test <- tibble(
    truth = fct_relevel(test$ar_primary, "0","1"),
    cal_score = calibrate_preds(fit_final, cal_fit, test)
  ) %>%
    mutate(pred_cls = factor(ifelse(cal_score >= thr, "1","0"), levels=c("0","1")))
  
  tibble(
    model = model_name,
    seed = seed,
    acc  = yardstick::accuracy(preds_test, truth, pred_cls)$.estimate,
    ba   = yardstick::bal_accuracy(preds_test, truth, pred_cls)$.estimate,
    recall = yardstick::sens(preds_test, truth, pred_cls, event_level = "second")$.estimate,
    spec   = yardstick::spec(preds_test, truth, pred_cls, event_level = "second")$.estimate,
    precision = yardstick::precision(preds_test, truth, pred_cls, event_level = "second")$.estimate,
    f1   = yardstick::f_meas(preds_test, truth, pred_cls, event_level = "second")$.estimate,
    rocA = yardstick::roc_auc_vec(truth = preds_test$truth, estimate = preds_test$cal_score, event_level="second"),
    prA  = yardstick::pr_auc_vec (truth = preds_test$truth, estimate = preds_test$cal_score, event_level="second")
  )
}


##################### Multi-seed Run #####################
random_seeds <- c(44542169,162339537,296891153,2080800158,
                  825781625,6315874,1525767726,1579953055,
                  1794694103,1214191569)

rf_results <- map_dfr(random_seeds, ~ evaluate_model_with_seed(rf_tune, "Random Forest", .x, data, rec))
xgb_results <- map_dfr(random_seeds, ~ evaluate_model_with_seed(boost_tune, "XGBoost", .x, data, rec))


all_results <- bind_rows(rf_results, xgb_results)

summary_results <- all_results %>%
  group_by(model) %>%
  summarise(across(c(acc, ba, recall, spec, precision, f1, rocA, prA),
                   list(mean = mean, median = median, sd = sd),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

print(all_results)
print(summary_results)

write.csv(summary_results, "result_rfxgb.csv")


##################### Identify Best Seed #####################
rf_summary <- rf_results %>%
  summarise(across(c(acc, ba, recall, f1, rocA, prA), mean, na.rm = TRUE))
print(rf_summary)

best_seed <- rf_results %>%
  arrange(desc(ba)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(seed)
cat("Best-performing seed:", best_seed, "\n")


######################## Retrain Best Model ########################
best_seed <- 162339537 
set.seed(best_seed)
spl <- initial_split(data, prop = 0.8, strata = "ar_primary")
train <- training(spl)
test  <- testing(spl)

rec_final <- prep(rec, training = train)

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(rf_tune)

params <- wf %>% extract_parameter_set_dials()
params <- params %>%
  update(mtry = finalize(mtry(), train %>% select(-ar_primary)))
folds <- vfold_cv(train, v = 5, strata = ar_primary)

res <- tune_bayes(
  wf, resamples = folds, param_info = params,
  initial = 10, iter = 20,
  metrics = metric_set(bal_accuracy),
  control = control_bayes(no_improve = 10, verbose = FALSE, parallel_over = "resamples")
)

best <- select_best(res, metric = "bal_accuracy")
final_wf <- finalize_workflow(wf, best)

fit_final_shap <- fit(final_wf, data = train)

# Extract model and recipe
rec_final     <- extract_recipe(fit_final_shap)
rf_model_shap <- extract_fit_engine(fit_final_shap)

rename_dict <- c(
  # --- Continuous / ratio variables ---
  "pir" = "Poverty–Income Ratio",
  "LBXEOPCT" = "Eosinophil % (blood)",
  "nlr" = "Neutrophil–Lymphocyte Ratio",
  "LBXIGE" = "Total IgE concentration",
  "crp_mg_l" = "C-Reactive Protein (mg/L)",
  "cot_ng_ml" = "Cotinine (ng/mL)",
  
  # --- Age groups ---
  "age_grp_X40.59" = "Age 40–59",
  "age_grp_X60.84" = "Age 60–84",
  
  # --- Demographics ---
  "sex_Female" = "Sex: Female",
  "race_Non.Hispanic.Black" = "Race: Non-Hispanic Black",
  "race_Non.Hispanic.White_x_LBXIGE" = "Race × IgE: Non-Hispanic White",
  
  # --- Home size ---
  "home_size_X2" = "Home size: 2 rooms",
  "home_size_X3" = "Home size: 3 rooms",
  "home_size_X4" = "Home size: 4 rooms",
  "home_size_X5" = "Home size: 5 rooms",
  "home_size_X7" = "Home size: 7+ rooms",
  
  # --- Education ---
  "edu_a_X9.11th..incl.12th.no.diploma." = "Education: 9–11th (no diploma)",
  "edu_a_High.school.GED" = "Education: High school/GED",
  "edu_a_Some.college.AA" = "Education: Some college/AA",
  "edu_a_College.graduate.or.above" = "Education: College graduate+",
  
  # --- Household income ---
  "house_income_X.10.000..14.999" = "Household income: $10–15k",
  "house_income_X.15.000..19.999" = "Household income: $15–20k",
  "house_income_X.20.000..24.999" = "Household income: $20–25k",
  "house_income_X.25.000..34.999" = "Household income: $25–35k",
  "house_income_X.35.000..44.999" = "Household income: $35–45k",
  "house_income_X.45.000..54.999" = "Household income: $45–55k",
  "house_income_X.55.000..64.999" = "Household income: $55–65k",
  "house_income_X.65.000..74.999" = "Household income: $65–75k",
  "house_income_X...75.000" = "Household income: ≥$75k",
  
  # --- Smoking status ---
  "smq_status_Former" = "Smoking: Former",
  "smq_status_Current" = "Smoking: Current",
  
  # --- BMI categories ---
  "bmi_cat_Overweight" = "BMI: Overweight",
  "bmi_cat_Obesity.I" = "BMI: Obesity I",
  "bmi_cat_Obesity.II" = "BMI: Obesity II",
  "bmi_cat_Obesity.III" = "BMI: Obesity III",
  "bmi_cat_unknown" = "BMI: Unknown",
  
  # --- Health and asthma ---
  "fam_asthma_ever_No" = "Family asthma history: No",
  "asthma_ever_No" = "Asthma history: No",
  "health_Very.Good" = "Health: Very good",
  "health_Good" = "Health: Good",
  "health_Fair" = "Health: Fair",
  
  # --- Allergies and pets ---
  "allergy_any_Yes" = "Any allergy: Yes",
  "animal_allergy_Yes" = "Animal allergy: Yes",
  "any_pet_Yes" = "Has pet: Yes",
  "mildew_Yes" = "Mildew exposure: Yes",
  
  # --- Home type ---
  "home_type_X1.family.house..detached" = "Home type: Detached house",
  "home_type_Apartment" = "Home type: Apartment",
  "home_type_Mobile.home.or.trailer" = "Home type: Mobile home/trailer",
  
  # --- Home year ---
  "home_year_X1950.to.1959" = "Home built: 1950–1959",
  "home_year_X1960.to.1977" = "Home built: 1960–1977",
  "home_year_X1978.to.1989" = "Home built: 1978–1989",
  "home_year_X1990.to.present" = "Home built: 1990–present",
  "home_year_Before.1940" = "Home built: Before 1940",
  "home_year_unknown" = "Home built: Unknown",
  
  # --- Number of rooms ---
  "num_rooms_X3" = "Rooms: 3",
  "num_rooms_X4" = "Rooms: 4",
  "num_rooms_X5" = "Rooms: 5",
  "num_rooms_X6" = "Rooms: 6",
  "num_rooms_X7" = "Rooms: 7",
  "num_rooms_X8" = "Rooms: 8",
  "num_rooms_X9" = "Rooms: 9",
  "num_rooms_other" = "Rooms: Other",
  
  # --- Insurance ---
  "any_ins_unknown" = "Insurance: Unknown",
  
  # --- Time of visit (seasonality proxy) ---
  "time_visit_X1" = "Visit time: Q1",
  "time_visit_X2.3" = "Visit time: Q2–Q3",
  "time_visit_X4.9" = "Visit time: Q4–Q9",
  "time_visit_X10.12" = "Visit time: Q10–Q12",
  "time_visit_X13." = "Visit time: Q13+",
  
  # --- Interactions ---
  "allergy_any_Yes_x_race_Non.Hispanic.White" = "Allergy × Race: Non-Hispanic White",
  "allergy_any_Yes_x_race_Non.Hispanic.Black" = "Allergy × Race: Non-Hispanic Black",
  "allergy_any_Yes_x_cot_ng_ml" = "Allergy × Cotinine"
)


###############################################################
train_baked <- bake(rec_final, new_data = train)
test_baked  <- bake(rec_final, new_data = test)

X_train_shap <- train_baked %>% select(-ar_primary)
X_test_shap  <- test_baked  %>% select(-ar_primary)

stopifnot(identical(names(X_test_shap),
                    rf_model_shap$forest$independent.variable.names))

###############################################################
##### SHAP computation
###############################################################
library(fastshap)
library(dplyr)
library(ggplot2)
library(ranger)
library(purrr)


random_seeds <- c(
  44542169,162339537,296891153,2080800158,
  825781625,6315874,1525767726,1579953055,
  1794694103,1214191569
)

for (s in random_seeds) {
  cat("▶️  Running SHAP for seed:", s, "...\n")
  set.seed(s)
  
  shap_values_shap_s <- fastshap::explain(
    object = rf_model_shap,
    X = X_test_shap,
    pred_wrapper = function(object, newdata)
     predict(object, data = newdata, type = "response")$predictions[, 2],
     nsim = 20
  )
  saveRDS(
    shap_values_shap_s,
    paste0("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/SHAP/shap_seed_", s, ".rds")
  )
  cat("✅  Finished seed:", s, "\n\n")
}

save_dir <- "/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/SHAP"
shap_files <- list.files(save_dir, pattern = "^shap_seed_.*\\.rds$", full.names = TRUE)
cat("📂  Detected", length(shap_files), "SHAP files.\n")

shap_list <- map(shap_files, readRDS)
dims <- unique(map(shap_list, dim))
cat("✅  Unique dimensions found:\n")
print(dims)

stopifnot(all(sapply(shap_list, function(x) all(dim(x) == dim(shap_list[[1]])))))
shap_values_shap_avg <- Reduce("+", shap_list) / length(shap_list)
avg_path <- file.path(save_dir, "shap_values_shap_avg.rds")
saveRDS(shap_values_shap_avg, avg_path)


################## Rename variable names for visualization
rename_columns <- function(df, dict) {
  colnames(df) <- ifelse(colnames(df) %in% names(dict),
                         dict[colnames(df)],
                         colnames(df))
  df
}

shap_values_shap_avg<-readRDS("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/SHAP/shap_values_shap_avg.rds")
shap_values_viz <- rename_columns(as.data.frame(shap_values_shap_avg), rename_dict)
test_baked_viz  <- rename_columns(as.data.frame(X_test_shap), rename_dict)

stopifnot(identical(colnames(shap_values_viz), colnames(test_baked_viz)))
cat("✅ Renamed successfully for visualization.\n")

shap_values_mat <- as.matrix(shap_values_viz)
saveRDS(test_baked, "test_baked.rds")
saveRDS(rename_dict, "rename_dict.rds")


# ============================================================
# 📦 Save essential SHAP and model objects for Force Dashboard
# ============================================================

# Define output directory
out_dir <- "/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/programs/support/force_dash1"

# Create the folder if it does not exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save key model and SHAP objects
saveRDS(rf_model_shap, file = file.path(out_dir, "rf_model_shap.rds"))
saveRDS(shap_values_mat, file = file.path(out_dir, "shap_values_mat.rds"))
saveRDS(shap_values_viz, file = file.path(out_dir, "shap_values_viz.rds"))
saveRDS(shap_values_shap_avg, file = file.path(out_dir, "shap_values_shap_avg.rds"))
saveRDS(test_baked_viz, file = file.path(out_dir, "test_baked_viz.rds"))
saveRDS(test, file = file.path(out_dir, "test_data.rds"))
saveRDS(X_test_shap, file = file.path(out_dir, "X_test_shap.rds"))

cat("✅ All key SHAP and model objects successfully saved to:", out_dir, "\n")


#####################################################################
out_dir2 <- "/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/programs/support/force_dash2"

saveRDS(shap_values_shap_avg, file = file.path(out_dir2, "shap_values_shap_avg.rds"))
saveRDS(test, file = file.path(out_dir2, "test.rds"))
saveRDS(test_baked, file = file.path(out_dir2, "test_baked.rds"))





