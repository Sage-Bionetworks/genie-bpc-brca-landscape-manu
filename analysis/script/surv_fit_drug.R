# Description: Fit the predictors of survival from distant metastasis.
# Additionally, fit these for hormone receptor subtypes.

n_boot <- 300 
boot_draw_seed <- 102039

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls('R'), .f = source)

read_folder <- here('data', 'survival', 'drug', 'prepared_data')
# create the folder if it does not exist:
fs::dir_create(here('data', 'survival', 'drug', 'fit_outputs'))
write_folder <- here('data', 'survival', 'drug', 'fit_outputs')

future::plan(strategy = multisession, workers = 6) 

tic()

dft_drug <- readr::read_rds(
  file = here(read_folder, 'dmet_surv_cdk.rds')
)

# in order to include bca_subtype_simple_f as a covariate we will need dummies 
#  for it:
dft_drug <- dft_drug %>%
  # make a simpler, standard column name version:
  mutate(
    bca = case_when(
      bca_subtype_f_simple %in% "HR+, HER2-" ~ "hr_pos_her2_neg",
      bca_subtype_f_simple %in% "HER2+" ~ "her2_pos",
      bca_subtype_f_simple %in% "Triple Negative" ~ "trip_neg",
      bca_subtype_f_simple %in% "(NC or NR)" ~ "nc_nr",
      T ~ "error"
    )
  ) %>%
  dummy_cols(., select_columns = "bca") %>% 
  select(-bca, -bca_hr_pos_her2_neg) 




# These are the default adjustments - but I like being explicit.
vec_dmet_surv_confounders <- c("age_dx_c", 
                               "stage_dx_iv_num",
                               "birth_year_c",
                               "white",
                               "hispanic")

boot_models_cdk <- surv_fit_dmet_wrap_no_lt(
  dat = dft_drug,
  boot_rep = n_boot, 
  main_seed = boot_draw_seed,
  additional_features = character(0),
  outcome_var = "tt_os_drug_yrs",
  outcome_ind = "os_dx_status"
) 



# Save all the fitted models as RDS:
readr::write_rds(
  x = boot_models_cdk,
  file = here(write_folder, 'fit_dmet_cdk.rds')
)

