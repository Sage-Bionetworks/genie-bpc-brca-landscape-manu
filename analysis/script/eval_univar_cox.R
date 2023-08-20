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



eval_beta_bias <- function(true_beta_valid, coef_dat) {
  beta_nonzero <- true_beta_valid[abs(true_beta_valid) > 0.001]
  
  coef_dat %<>%
    filter(term %in% names(beta_nonzero)) %>%
    select(term, estimate)
  
  beta_dat <- tibble(
    term = names(beta_nonzero),
    truth = beta_nonzero
  )
  
  bias_dat <- left_join(
    beta_dat,
    coef_dat,
    by = "term",
    relationship = "one-to-one"
  )
  
  if (any(is.na(bias_dat$estimate))) {
    cli::cli_alert_warning("Some terms were not found in coef_dat")
  }
  
  bias_dat %<>% mutate(diff = estimate - truth)
  
  bias_dat %<>%
    summarize(
      avg_abs_bias = mean(abs(diff), na.rm = T),
      avg_bias = mean(diff, na.rm = T),
      se_bias = sd(diff, na.rm = T)
    )
  
  if (nrow(bias_dat) > 1) {
    cli::cli_abort("Something went wrong, return dataframe has multiple rows.")
  }
  
  return(bias_dat) # returns a dataframe, so you'll have to unpack it.
  
}

eval_beta_bias(true_beta_valid = test_beta,
               coef_dat = test_coef_dat)

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



# next steps: 
# - calculating the average bias.  
# - Probably calculating versions of the beta and generated data coefficients that are filtered for at least one positive would be good too.
# - adding straight ridge and lasso.  Ridge won't have AUC and that's fine.

