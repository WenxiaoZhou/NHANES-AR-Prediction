set.seed(20250928)
load("/Users/zhouwenxiao/Desktop/ANLY699/Final Project Related/data/data.RData")
data <- data_model

library(glmnet)

# (2) Split ----------
spl   <- initial_split(data, prop = 0.8, strata = ar_primary)
train <- training(spl)

# generate the interaction terms
# (.)^2 = all second-order interactions
train_complete <- train %>% drop_na()
X <- model.matrix(ar_primary ~ (.)^2, data = train_complete)[,-1]
y <- train_complete$ar_primary

# LASSO logistic regression
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1)
# check what interaction terms are kept
coefs <- coef(cvfit, s = "lambda.min")
active_vars <- rownames(coefs)[which(coefs != 0)]
active_coefs <- coefs[which(coefs != 0), ]
active_interactions <- tibble(variable = active_vars, coef = as.numeric(active_coefs)) %>%
  filter(str_detect(variable, ":"))

