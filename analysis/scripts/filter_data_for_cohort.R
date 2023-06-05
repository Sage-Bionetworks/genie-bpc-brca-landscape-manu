# Description:  Filters the data to only index cancer rows.  Also removes 
#  participants known to have Sarcoma, a cohort criterion agreed upon by group.

library(dplyr)
library(purrr)
library(readr)
library(magrittr)

# Load all helper functions
purrr::walk(.x = fs::dir_ls('R'), .f = source)

dft_pt <- readr::read_csv(here('data-raw', 'patient_level_dataset.csv'))
dft_ca_ind <- readr::read_csv(here('data-raw', 'cancer_level_dataset_index.csv'))
dft_cpt <- readr::read_csv(here('data-raw', 'cancer_panel_test_level_dataset.csv'))
dft_img <- readr::read_csv(here('data-raw', 'imaging_level_dataset.csv'))
dft_med_onc <- readr::read_csv(here('data-raw', 'med_onc_note_level_dataset.csv'))
dft_path <- readr::read_csv(here('data-raw', 'pathology_report_level_dataset.csv'))
dft_reg <- readr::read_csv(here('data-raw', 'regimen_cancer_level_dataset.csv'))
dft_tm <- readr::read_csv(here('data-raw', 'tm_level_dataset.csv'))

# Not touching some datasets including:
# - cancer_level_dataset_non_index.csv
# - all genetic datasets and gene panel level datasets (processed separately)

# Note:  ca_hist_adeno_squamous is about 20% missing, so this criterion has
# limited appeal.
vec_not_sarcoma <- dft_ca_ind %>%
  select(record_id, ca_seq, ca_hist_adeno_squamous) %>%
  filter(!(ca_hist_adeno_squamous %in% "Sarcoma")) %>%
  pull(record_id)

# Some of these datasets can be filtered by this vector:
dft_pt %<>% filter(record_id %in% vec_not_sarcoma)
dft_ca_ind %<>% filter(record_id %in% vec_not_sarcoma)
dft_img %<>% filter(record_id %in% vec_not_sarcoma)
dft_med_onc %<>% filter(record_id %in% vec_not_sarcoma)
dft_path %<>% filter(record_id %in% vec_not_sarcoma)
dft_tm %<>% filter(record_id %in% vec_not_sarcoma)

dft_cohort_keys <- dft_ca_ind %>% select(record_id, ca_seq)
chk_keys_unique <- count(dft_cohort_keys, record_id, ca_seq) %>%
  pull(n) %>% 
  max %>%
  is_in(1)
if (!chk_keys_unique) {
  cli::cli_abort("Duplicate keys found in filter_data_for_cohort.R")
}

key_filt_help <- function(dat) {
  left_join(
    dft_cohort_keys,
    dat,
    by = c("record_id", "ca_seq"),
    multiple = "all" # default behavior in SQL and dplyr - just silences.
  )
}

dft_cpt %<>% key_filt_help(.)
dft_reg %<>% key_filt_help(.)

fs::dir_create(here('data', 'clin_data_cohort'))

readr::write_rds(
  x = dft_pt,
  file = here('data', 'clin_data_cohort', 'dft_pt.rds')
)


# Not a best practice to use the names, but it will work.
write_help <- function(dat) {
  nm <- deparse(substitute(dat))
  readr::write_rds(
    x = dat,
    file = here('data', 'clin_data_cohort', 
                paste0(nm,'.rds'))
  )
}

write_help(dft_pt)
write_help(dft_ca_ind)
write_help(dft_img)
write_help(dft_med_onc)
write_help(dft_path)
write_help(dft_tm)
write_help(dft_cpt)
write_help(dft_reg)




