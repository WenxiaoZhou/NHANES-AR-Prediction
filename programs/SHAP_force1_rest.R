
####### Group-level SHAP aggregation
# Mean SHAP Heatmap for Allergy_any: Yes x Race x Cotinine Interaction
library(dplyr)
library(ggplot2)
group_shap <- as.data.frame(shap_values_shap_avg) %>%
  mutate(
    race = test$race,  # categorical race variable
    cotinine_group = cut(test$cot_ng_ml, breaks = c(-Inf, 1, 10, Inf))  # exposure levels
  ) %>%
  group_by(race, cotinine_group) %>%
  summarise(
    mean_shap_allergy = mean(allergy_any_Yes_x_cot_ng_ml, na.rm = TRUE),
    .groups = "drop"
  )
ggplot(group_shap %>% filter(!is.na(cotinine_group)),
       aes(x = cotinine_group, y = race, fill = mean_shap_allergy)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#3182bd",    # blue: decreases predicted probability
    mid = "white",
    high = "#e34a33",   # red: increases predicted probability
    midpoint = 0,
    name = "Mean SHAP\n(Allergy = Yes)"
  ) +
  labs(
    title = "Mean SHAP for 'Allergy Any × Cotinine Group' across Race",
    x = "Cotinine exposure group (ng/mL)",
    y = "Race"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    panel.grid = element_blank()
  )

#sv_dependence(sv, v = "allergy_any_Yes", color_var = "cot_ng_ml")

# ============================================================
# Interpretation guide:
#   🔵 Blue  = allergy history decreases predicted AR probability
#   🔴 Red   = allergy history increases predicted AR probability
#   ⚪ White = near-neutral SHAP contribution
# ============================================================

# Mean SHAP Heatmap for Race x IgE Interaction
df_IgE <- as.data.frame(shap_values_shap_avg) %>%
  mutate(
    race = test$race,
    LBXIGE = test$LBXIGE,
    IgE_group = cut(LBXIGE, breaks = c(-Inf, 1, 3, 10, Inf),
                    labels = c("<=1", "1–3", "3–10", ">10"),
                    right=TRUE,
                    include.lowest = TRUE)
  )

# Compute mean SHAP by race × IgE group
group_IgE <- df_IgE %>%
  filter(!is.na(IgE_group)) %>%
  group_by(race, IgE_group) %>%
  summarise(mean_shap = mean(race_Non.Hispanic.White_x_LBXIGE, na.rm = TRUE),
            .groups = "drop")

# Detect sign reversal
if (range(group_IgE$mean_shap, na.rm = TRUE)[1] < 0 &
    range(group_IgE$mean_shap, na.rm = TRUE)[2] > 0) {
  cat("⚠️ SHAP direction reverses across race × IgE groups\n")
} else {
  cat("✅ SHAP direction consistent across groups\n")
}

ggplot(group_IgE, aes(IgE_group, race, fill = mean_shap)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#e34a33", midpoint = 0) +
  labs(
    title = "Mean SHAP for Race × IgE Interaction",
    x = "IgE group (ng/mL)", y = "Race", fill = "Mean SHAP"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12)
  )



##############################################################
##############################################################
##############################################################

#summary plot
library(ranger) 
library(ggplot2) 
library(plotly) 
library(shapviz)
library(dplyr)
library(htmltools)

# ===============================================================
#  SHAPviz with renamed, human-readable feature names
# ===============================================================
library(ranger)
library(dplyr)
library(shapviz)
library(ggplot2)
library(readxl)
shap_values_mat <- as.matrix(shap_values_shap_avg)

rf_pred <- ranger:::predict.ranger(rf_model_shap, data = X_test_shap)
base_value <- mean(rf_pred$predictions[, 2])
cat("✅ Base value computed:", base_value, "\n")

test_baked_sv <- test_baked %>% select(-ar_primary)

# 2️⃣ Create SHAPviz object
sv <- shapviz(shap_values_mat, X = test_baked_sv, baseline = base_value)

# 3️⃣ Define manual rename map for key interactions
rename_map <- c(
  "allergy_any_Yes_x_cot_ng_ml" = "Any Allergy: Yes × Cotinine (ng/mL)",
  "race_Non.Hispanic.White_x_LBXIGE" = "Race: Non-Hispanic White × Total IgE",
  "allergy_any_Yes_x_race_Non.Hispanic.White" = "Any Allergy × Race: Non-Hispanic White",
  "allergy_any_Yes_x_race_Non.Hispanic.Black" = "Any Allergy × Race: Non-Hispanic Black"
)

# 4️⃣ Read dummy variable label list from Excel
dummy_labels <- read_excel("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/programs/support/force_dash2/dummy_namelist.xlsx")
colnames(dummy_labels)[1:4] <- c("dummy", "group", "display_old", "rename")

rename_df <- dummy_labels %>%
  select(dummy, rename) %>%
  mutate(rename = trimws(rename))

# 5️⃣ Build unified rename mapping
rename_combined <- rename_df$rename
names(rename_combined) <- rename_df$dummy
rename_combined <- c(rename_combined, rename_map)

# 6️⃣ Apply to both sv$X and sv$S to maintain matching column names
old_names <- colnames(sv$X)
new_names <- ifelse(
  old_names %in% names(rename_combined),
  rename_combined[old_names],
  old_names
)

colnames(sv$X) <- new_names
colnames(sv$S) <- new_names

# 7️⃣ Plot with readable names
sv_importance(sv, kind = "beeswarm", max_display = 20)

