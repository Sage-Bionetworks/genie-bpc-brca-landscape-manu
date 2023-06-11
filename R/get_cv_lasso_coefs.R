get_cv_lasso_coefs <- function(las_fit) {
  las_fit %>%
    coef(., s = "lambda.min") %>%
    tidy_cv_glmnet(., exp_coef = T, remove_zero = F)
}