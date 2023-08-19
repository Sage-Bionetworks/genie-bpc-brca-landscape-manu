# Description: Runs the univariate cox method on simulated datasets.

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

sim_n80 <- readr::read_rds(
  here('sim', 'run_methods', 'gen_dat_one_n80_cox_uni.rds')
)

sim_n500 <- readr::read_rds(
  here('sim', 'run_methods', 'gen_dat_one_n500_cox_uni.rds')
)

# Example of running one time:
test_beta <- sim_n80 %>% slice(2) %>% pull(beta) %>% unlist(.)
test_coef_df <- sim_n80 %>% slice(2) %>% pull(coef_est) %>% `[[`(.,1)
test_gen_dat <- sim_n80 %>% slice(2) %>% pull(gen_dat) %>% `[[`(.,1)
# 
# split_decision_values(
#   true_beta = test_beta,
#   coef_dat = test_coef_df,
#   d_name = "p.value"
# )
eval_beta_auc_pval(
  true_beta = test_beta,
  coef_dat = test_coef_df
)


eval_beta_bias <- function(true_beta, coef_dat) {
  true_beta <- 
  
  
}

# Add AUC values in:
sim_n80 %<>%
  mutate(
    analysis_method = "method_univar_cox",
    auc = purrr::map2_dbl(
      .x = beta,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_auc_pval(
          true_beta = b,
          coef_dat = c,
          return_type = "auc"
        )
      })
    )
  )

sim_n500 %<>%
  mutate(
    analysis_method = "method_univar_cox",
    auc = purrr::map2_dbl(
      .x = beta,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_auc_pval(
          true_beta = b,
          coef_dat = c,
          return_type = "auc"
        )
      })
    )
  )

# next steps: 
# - calculating the average bias.  
# - Probably calculating versions of the beta and generated data coefficients that are filtered for at least one positive would be good too.
# - adding straight ridge and lasso.  Ridge won't have AUC and that's fine.

