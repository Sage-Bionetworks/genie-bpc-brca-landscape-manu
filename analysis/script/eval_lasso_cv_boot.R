# Description: Runs the univariate cox method on simulated datasets.

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

cli::cli_alert_danger("Loading the first 200 reps only (pilot)")

sim_n80 <- readr::read_rds(
  here('sim', 'run_methods', 'gen_dat_one_n80_lasso_cv_boot_f200.rds')
)

sim_n500 <- readr::read_rds(
  here('sim', 'run_methods', 'gen_dat_one_n500_lasso_cv_boot_f200.rds')
)


# sim_n80 %>% slice(1) %>% pull(coef_est) %>% `[[`(.,1)


# There's an if statement here because eventually I'd like this to be done at the previous step, but I didn't have time to add it previously.
# For now it just takes a few seconds to compute it here, so off we go:
if (!("coef_est" %in% colnames(sim_n80))) {
  sim_n80 %<>%
    mutate(
      coef_est = purrr::map(
        .x = fit,
        .f = (function(f) resample_coef_helper(f))
      )
    )
  
  sim_n500 %<>%
    mutate(
      coef_est = purrr::map(
        .x = fit,
        .f = (function(f) resample_coef_helper(f))
      )
    ) 
  
}



###########################
# Add in AUC calculation: #
###########################

sim_n80 %<>%
  mutate(
    auc = purrr::map2_dbl(
      .x = beta_valid,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_auc_stability(
          true_beta = b,
          coef_dat = c,
          return_type = "auc"
        )
      })
    )
  )

sim_n500 %<>%
  mutate(
    auc = purrr::map2_dbl(
      .x = beta_valid,
      .y = coef_est,
      .f = (function(b, c) {
        eval_beta_auc_stability(
          true_beta = b,
          coef_dat = c,
          return_type = "auc"
        )
      })
    )
  )



# test_beta <- sim_n80 %>% slice(1) %>% pull(beta_valid) %>% `[[`(.,1) 
# test_coef <- sim_n80 %>% slice(1) %>% pull(coef_est) %>% `[[`(.,1) 
# 
# setdiff(names(test_beta), test_coef$term)
# setdiff(test_coef$term, names(test_beta))


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


# A good test case on bias:
# test_beta <- sim_n500 %>% slice(1) %>% pull(beta_valid) %>% `[[`(.,1)
# test_beta <- test_beta[abs(test_beta) > 0.0000001] 
# test_coef_est <- sim_n500 %>% slice(1) %>% pull(coef_est) %>% `[[`(.,1)
# test_coef_est %>%
#   filter(term %in% names(test_beta)) %>%
#   mutate(truth = test_beta,
#          diff = estimate - truth) %>%
#   summarize(
#     abs_bias = mean(abs(diff)),
#     bias = mean(diff)
#   )
# # should match:
# sim_n500 %>% slice(1) %>% select(contains("bias"))
