# Description: Runs the univariate cox method on simulated datasets.

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

sim_n80 <- readr::read_rds(
  here('sim', 'run_methods',
       'gen_dat_one_n80_lasso_5fcv.rds')
)

sim_n500 <- readr::read_rds(
  here('sim', 'run_methods', 
       'gen_dat_one_n500_lasso_5fcv.rds')
)

sim_n80 %<>% unnest(fit)
sim_n500 %<>% unnest(fit)



###########################
# Add in AUC calculation: #
###########################
# 
# sim_n80 %<>%
#   mutate(
#     auc = purrr::map2_dbl(
#       .x = beta_valid,
#       .y = coef_est,
#       .f = (function(b, c) {
#         eval_beta_auc_stability(
#           true_beta = b,
#           coef_dat = c,
#           return_type = "auc"
#         )
#       })
#     )
#   )
# 
# sim_n500 %<>%
#   mutate(
#     auc = purrr::map2_dbl(
#       .x = beta_beta,
#       .y = coef_est,
#       .f = (function(b, c) {
#         eval_beta_auc_stability(
#           true_beta = b,
#           coef_dat = c,
#           return_type = "auc"
#         )
#       })
#     )
#   )




############################
# Add in bias calculation: #
############################
sim_n80 %<>%
  # add several bias metrics in:
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
  # add several bias metrics in:
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
  file = here('sim', 'evaled_methods', 'gen_dat_one_n80_lasso_cv_boot_f200.rds')
)

readr::write_rds(
  x = sim_n500,
  file = here('sim', 'evaled_methods', 'gen_dat_one_n500_lasso_cv_boot_f200.rds')
)


