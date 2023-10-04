# Description: Process the raw model fits into results.
# Author: Alex Paynter



library(fs); library(purrr); library(here);
purrr::walk(.x = fs::dir_ls(here("R")), .f = source)


raw_folder <- here("data", "survival", 'drug', 'fit_outputs')
# Create the folder if it does not exist:
fs::dir_create('data', 'survival', 'drug', 'fit_outputs', 'fit_summary')
output_folder <- here('data', 'survival', 'drug', 'fit_outputs', 'fit_summary')


boot_models_dmet_cdk <- readr::read_rds(
  here(raw_folder, "fit_dmet_cdk.rds")
)

resample_dmet_wrap <- function(dat) {
  resample_coef_helper(
    dat = dat,
    exp_coefs = F,
    estimate_col_rename = "log_hr" 
  )
}

dft_coef_cdk <- resample_dmet_wrap(
  dat = boot_models_dmet_cdk
)


# Write the outputs:
write_wrap_surv_sum <- function(obj, name) {
    readr::write_rds(
      x = obj,
      file = here(output_folder, paste0(name, ".rds"))
    )
}

write_wrap_surv_sum(
  dft_coef_cdk,
  "coef_cdk"
)


