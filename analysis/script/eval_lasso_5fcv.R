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

# test_beta <- sim_n500 %>% slice(1) %>% pull(beta_valid) %>% `[[`(.,1)
# test_dat <- sim_n500 %>% slice(1) %>% pull(dropout_dat) %>% `[[`(.,1)
# eval_beta_auc_lambda(test_beta, test_dat, return_type = "plot")

###########################
# Add in AUC calculation: #
###########################

sim_n80 %<>%
  mutate(
    auc = purrr::map2_dbl(
      .x = beta_valid,
      .y = dropout_dat,
      .f = (function(b, d) {
        eval_beta_auc_lambda(
          true_beta = b,
          dropout_dat = d,
          return_type = "auc"
        )
      })
    )
  )

sim_n500 %<>%
  mutate(
    auc = purrr::map2_dbl(
      .x = beta_valid,
      .y = dropout_dat,
      .f = (function(b, d) {
        eval_beta_auc_lambda(
          true_beta = b,
          dropout_dat = d,
          return_type = "auc"
        )
      })
    )
  )


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
        # I didn't save the names right on the first run:
        if ("log_hr" %in% names(c)) {
          c %<>% rename(estimate = log_hr)
        }
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
        # I didn't save the names right on the first run:
        if ("log_hr" %in% names(c)) {
          c %<>% rename(estimate = log_hr)
        }
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
  file = here('sim', 'evaled_methods', 'gen_dat_one_n80_lasso_5fcv.rds')
)

readr::write_rds(
  x = sim_n500,
  file = here('sim', 'evaled_methods', 'gen_dat_one_n500_lasso_5fcv.rds')
)


