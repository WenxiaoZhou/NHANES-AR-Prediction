library(tidyverse)
library(reshape2)
library(plotly)
library(ggthemes)

# Compute SHAP summary per feature
shap_values_shap <- as.data.frame(shap_values_shap)
shap_summary <- shap_values_shap %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "shap") %>%
  group_by(feature) %>%
  summarise(
    mean_shap = mean(shap, na.rm = TRUE),
    mean_abs_shap = mean(abs(shap), na.rm = TRUE),
    pos_ratio = mean(shap > 0, na.rm = TRUE),
    neg_ratio = mean(shap < 0, na.rm = TRUE)
  ) %>%
  mutate(
    reversal_index = pmin(pos_ratio, neg_ratio) / (pos_ratio + neg_ratio),
    feature_clean = ifelse(feature %in% names(rename_dict),
                           rename_dict[feature], feature)
  ) %>%
  arrange(desc(reversal_index))

# Top 5 features
top5_reversal <- shap_summary %>% slice_max(reversal_index, n = 5)
knitr::kable(top5_reversal, digits = 3, caption = "Top 5 Features with Highest SHAP Direction Reversal")
