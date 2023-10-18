# Description: Process the raw model fits into results.
# Author: Alex Paynter

library(fs); library(purrr); library(here);
purrr::walk(.x = fs::dir_ls(here("R")), .f = source)

raw_folder <- here("data", "survival", 'drug', 'fit_outputs')
# now that I'm using nested models this is totally tractable:
output_folder <- here('data', 'survival', 'drug', 'fit_outputs')

dft_drug_surv_nest <- readr::read_rds(
  here(raw_folder, "fit_dmet_drug_nest.rds")
)

resample_dmet_wrap <- function(dat) {
  resample_coef_helper(
    dat = dat,
    exp_coefs = F,
    estimate_col_rename = "log_hr" 
  )
}

dft_drug_surv_nest %<>%
  mutate(
    coef_os = purrr::map(
      .x = boots_os,
      .f = resample_dmet_wrap
    ),
    coef_os_no_conf = purrr::map(
      .x = boots_os_no_conf,
      .f = resample_dmet_wrap
    ),
    
    coef_pfs = purrr::map(
      .x = boots_pfs,
      .f = resample_dmet_wrap
    ),
    coef_pfs_no_conf = purrr::map(
      .x = boots_pfs_no_conf,
      .f = resample_dmet_wrap
    )
  )



# dft_drug_surv_nest %>% slice(1) %>% pull(coef_os) %>% `[[`(.,1) %>%
#   arrange(desc(stability))

# Glorious nested data makes me feel OK imbedding the plots right here:




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


