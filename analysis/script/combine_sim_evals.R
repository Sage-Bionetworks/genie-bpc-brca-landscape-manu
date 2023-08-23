# Description: Combines the simulation runs into dataframes for display and analysis.

library(fs); library(here); library(purrr);
purrr::walk(.x = fs::dir_ls('R'), .f = source)


read_wrap <- function(f) {
  readr::read_rds(
    here(
      'sim', 'evaled_methods', f
    )
  )
}

n80_cox <- read_wrap('gen_dat_one_n80_cox_uni.rds')
n500_cox <- read_wrap('gen_dat_one_n500_cox_uni.rds')

n80_lasso_5fcv <- read_wrap('gen_dat_one_n80_lasso_5fcv.rds')
n500_lasso_5fcv <- read_wrap('gen_dat_one_n500_lasso_5fcv.rds')

n80_lasso_cv_boot <- read_wrap(
  'gen_dat_one_n80_lasso_cv_boot_f200.rds' 
)
n500_lasso_cv_boot <- read_wrap(
  'gen_dat_one_n80_lasso_cv_boot_f200.rds' 
)


n80_cox %<>% mutate(n = 80)
n500_cox %<>% mutate(n = 500)

n80_lasso_5fcv %<>% mutate(n = 80)
n500_lasso_5fcv %<>% mutate(n = 500)

n80_lasso_cv_boot %<>% mutate(n = 80)
n500_lasso_cv_boot %<>% mutate(n = 500)

cli::cli_alert_danger(
  text = "Temp fix:  Cutting all models down to the first 200 runs for matching"
)
n80_cox %<>% slice(1:200)
n500_cox %<>% slice(1:200)
n80_lasso_5fcv %<>% slice(1:200)
n500_lasso_5fcv %<>% slice(1:200)


sim_sum_all <- bind_rows(
  n80_cox,
  n500_cox,
  n80_lasso_cv_boot,
  n500_lasso_cv_boot
)