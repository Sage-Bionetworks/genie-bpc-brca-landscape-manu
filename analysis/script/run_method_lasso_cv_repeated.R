# Description: Runs the method on the datasets

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

sim_n80 <- readr::read_rds(
  here('sim', 'gen_data', 'gen_dat_one_n80.rds')
)

sim_n500 <- readr::read_rds(
  here('sim', 'gen_data', 'gen_dat_one_n500.rds')
)

tic()

# Example of running one dataset:
test_dat <- sim_n500 %>%
  slice(1) %>%
  pull(gen_dat) %>%
  `[[`(1)

method_lasso_cv_repeat(
  test_dat,
  main_seed = 2019,
  rep = 100 
)

toc()
# # 
# # method_univar_cox(test_dat)
# 
# 
# sim_n80 %<>%
#   mutate(
#     analysis_method = "method_univar_cox",
#     gen_dat = purrr::map(
#       .x = gen_dat,
#       .f = (function(d) {
#         method_univar_cox(d, ignore_cols = "id_obs")
#       })
#     )
#   )
# 
# sim_n500 %<>%
#   mutate(
#     analysis_method = "method_univar_cox",
#     gen_dat = purrr::map(
#       .x = gen_dat,
#       .f = (function(d) {
#         method_univar_cox(d, ignore_cols = "id_obs")
#       })
#     )
#   )


# cli::cli_alert_success(
#   glue(
#     "Ran the cox univariate method in {toc()$callback_msg}."
#   )
# )
# 
# readr::write_rds(
#   x = sim_n80,
#   file = here('sim', 'run_methods', 'gen_dat_one_n80_lasso_cv_repeated.rds')
# )
# 
# readr::write_rds(
#   x = sim_n500,
#   file = here('sim', 'run_methods', 'gen_dat_one_n500_lasso_cv_repeated.rds')
# )


  
    