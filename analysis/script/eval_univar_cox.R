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
test_beta <- sim_n80 %>% slice(2) %>% pull(beta_valid) %>% unlist(.)
test_coef_dat <- sim_n80 %>% slice(2) %>% pull(coef_est) %>% `[[`(.,1)
test_gen_dat <- sim_n80 %>% slice(2) %>% pull(gen_dat_valid) %>% `[[`(.,1)
# 
# split_decision_values(
#   true_beta = test_beta,
#   coef_dat = test_coef_dat,
#   d_name = "p.value"
# )
# eval_beta_auc_pval(
#   true_beta = test_beta,
#   coef_dat = test_coef_dat
# )
# eval_beta_bias(true_beta_valid = test_beta,
#                coef_dat = test_coef_dat)

# Add AUC values in:
sim_n80 %<>%
  mutate(
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

# add several bias metrics in:
sim_n80 %<>%
  mutate(
    bias_dat = purrr::map2(
      .x = beta_valid,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_bias(
          true_beta_valid = b,
          coef_dat = c
        )
      })
    )
  ) %>%
  unnest(bias_dat)



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

# add several bias metrics in:
sim_n500 %<>%
  mutate(
    bias_dat = purrr::map2(
      .x = beta_valid,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_bias(
          true_beta_valid = b,
          coef_dat = c
        )
      })
    )
  ) %>%
  unnest(bias_dat)


readr::write_rds(
  x = sim_n80,
  file = here('sim', 'evaled_methods', 'gen_dat_one_n80_cox_uni.rds')
)


readr::write_rds(
  x = sim_n500,
  file = here('sim', 'evaled_methods', 'gen_dat_one_n500_cox_uni.rds')
)



