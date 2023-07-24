# Description: Fit the predictors of survival from distant metastasis.
# Additionally, fit these for hormone receptor subtypes.

n_boot <- 300
boot_draw_seed <- 102039

library(magrittr)
library(dplyr)
library(janitor)
library(here)
library(purrr)
library(fs)
library(ggplot2)
library(tidyr)
library(stringr)
library(survival)
library(glmnet)
library(tictoc)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

tic()

dft_dmet_surv <- readr::read_rds(
  file = here('data', 'survival', 'prepared_data', 'surv_dmet.rds')
)

# We will repeat this for four cohorts (including dft_dmet_surv itself):
dft_dmet_surv_trip_neg <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "Triple Negative")
dft_dmet_surv_hr_pos_her2_neg <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "HR+, HER2-")
dft_dmet_surv_her2_pos <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "HER2+")


# Big wrapper that creates the surv object, fits the survival, returns the datasets.
# Because the model is really specific to this one use case, we didn't put this
# wrapper into the /R directory.
surv_fit_dmet_wrap <- function(dat, boot_rep, main_seed) {
  y_dmet_os <- with(
    dat,
    Surv(
      time = tt_cpt_dmet_yrs_pos,
      time2 = tt_os_dmet_yrs,
      event = os_dx_status
    )
  )
  
  x_dmet_os <- dat %>%
    # just to get the output in readable order we do an alphabet sort first:
    select(order(colnames(dat))) %>%
    select(
      matches("_mut$"),
      matches("_cna$"),
      matches("_fus$"),
      age_dx_c,
      stage_dx_iv_num,
      birth_year_c,
      white,
      hispanic
    ) %>%
    as.matrix(.)
  
  set.seed(main_seed) # for drawing each boot seed below
  dft_boot_dmet_os <- tibble(boot_ind = 1:boot_rep)%>%
    mutate(boot_seed = sample.int(n = 10^7, size = n()))
  
  rtn <- dft_boot_dmet_os %>%
    mutate(
      fit = purrr::map(
        .x = boot_seed,
        .f = (function(x) {
          fit_boot_cv_ltrc_lasso(
            x_mat = x_dmet_os,
            y_df = dft_dmet_surv,
            y_t = "tt_cpt_dmet_yrs_pos",
            y_t2 = "tt_os_dmet_yrs",
            y_e = "os_dx_status",
            boot_seed = x,
          )
        })
      )
    )
  
  return(rtn)
}

boot_models_dmet_all <- surv_fit_dmet_wrap(
  dat = dft_dmet_surv, 
  boot_rep = n_boot, 
  main_seed = boot_draw_seed
) 

boot_models_dmet_trip_neg <- surv_fit_dmet_wrap(
  dat = dft_dmet_surv_trip_neg,
  boot_rep = n_boot, 
  main_seed = boot_draw_seed
) 

boot_models_dmet_hr_pos_her2_neg <- surv_fit_dmet_wrap(
  dat = dft_dmet_surv_hr_pos_her2_neg,
  boot_rep = n_boot, 
  main_seed = boot_draw_seed
) 

boot_models_dmet_her2_pos <- surv_fit_dmet_wrap(
  dat = dft_dmet_surv_her2_pos,
  boot_rep = n_boot, 
  main_seed = boot_draw_seed + 1
) 


vec_time_elapsed <- toc()
cli::cli_progress_message(
  paste0(
   vec_time_elapsed$callback_msg,
    " to run dmet boots with ",
    n_boots,
    " repetitions."
  )
)

# Save all the fitted models as RDS:
readr::write_rds(
  x = boot_models_dmet_all,
  file = here('data', 'survival', 'fit_outputs', 'fit_dmet_all.rds')
)
readr::write_rds(
  x = boot_models_dmet_trip_neg,
  file = here('data', 'survival', 'fit_outputs', 'fit_dmet_trip_neg.rds')
)
readr::write_rds(
  x = boot_models_dmet_hr_pos_her2_neg,
  file = here('data', 'survival', 'fit_outputs', 'fit_dmet_hr_pos_her2_neg.rds')
)
readr::write_rds(
  x = boot_models_dmet_her2_pos,
  file = here('data', 'survival', 'fit_outputs', 'fit_dmet_her2_pos.rds')
)

