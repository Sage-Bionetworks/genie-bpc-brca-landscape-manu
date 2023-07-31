# Description: Create the survival dataset indexing from the time of 
#   distant metastasis.

library(magrittr)
library(dplyr)
library(janitor)
library(yaml)
library(here)
library(purrr)
library(fs)
library(ggplot2)
library(tidyr)
library(stringr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)


# Load some data to play with:
dft_pt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_pt.rds')
)
dft_ca_ind <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
)
dft_cpt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_cpt.rds')
)
dft_gene_feat <- readr::read_rds(
  here('data', 'gene_feat_long.rds')
)

# Start of getting cohort info as we go:
dft_ca_ind %>%
  tabyl(bca_subtype_f_simple) %>%
  adorn_totals(.) %>%
  select(1:2) %>%
  mutate(
    dmet_filter = 0,
    entry_after_event_filter = 0
  )
         

dft_cpt %>%
  select(record_id, cpt_order_int, cpt_report_int) %>%
  filter(is.na(cpt_order_int) | is.na(cpt_report_int))
# Never missing - I've been lied to!



dft_dmet_timing <- get_dmet_timing(ca_ind_df = dft_ca_ind)

# Example:  Get all NGS which occur prior to metastatic cancer, or any which are
#   the first NGS test.
dft_cpt_dmet <- get_cpt_by_time(
  time_dat = dft_dmet_timing,
  time_var = "tt_y",
  cpt_dat = dft_cpt,
  always_keep_first = T
)

dft_cpt_dmet %<>%
  select(
    record_id, 
    ca_seq, 
    cpt_genie_sample_id,
    dx_cpt_rep_yrs,
    is_first_cpt,
    cpt_before_t,
    cpt_seq_assay_id
  )

# 
dft_clin_char <- dft_pt %>%
  mutate(
    white = case_when(
      is.na(naaccr_ethnicity_code) ~ NA_real_,
      naaccr_race_code_primary %in% "White" ~ 1,
      T ~ 0
    ),
    hispanic = case_when(
      is.na(naaccr_ethnicity_code) ~ NA_real_,
      naaccr_ethnicity_code %in% "Non-Spanish; non-Hispanic" ~ 0,
      T ~ 1
    )
  ) %>%
  select(record_id, white, hispanic, birth_year)

dft_clin_char <- dft_ca_ind %>%
  select(
    record_id, 
    ca_seq, # works because we've already selected one row per person.
    contains("bca_subtype"), # several versions.
    age_dx, 
    stage_dx_iv,
    dmets_stage_i_iii,
    dx_to_dmets_yrs,
    os_dx_status, 
    tt_os_dx_yrs,
    pfs_i_and_m_adv_status,
    tt_pfs_i_and_m_adv_yrs
  ) %>%
  full_join(., dft_clin_char, by = "record_id") 

# A helper function which only has relevance to analyses from dmet:
add_specific_dmet_vars <- function(dat) {
  dat %<>% mutate(
    dat,
    tt_cpt_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ dx_cpt_rep_yrs - dx_to_dmets_yrs,
    ),
    # glmnet won't take negative times:
    tt_cpt_dmet_yrs_pos = if_else(tt_cpt_dmet_yrs < 0, 0, tt_cpt_dmet_yrs)
  )
  return(dat)
}




# Do the fir
dft_gene_comb_all <- combine_cpt_gene_feat(
  dat_gene_feat = dft_gene_feat_dmet,
  dat_cpt = dft_cpt_dmet
)

dft_dmet_surv_all <- combine_clin_gene(
  dat_gene = dft_gene_comb_all,
  dat_clin = dft_clin_char
) %>%
  add_specific_dmet_vars(.)

dft_dmet_surv_all %<>%
  filter(!(tt_cpt_dmet_yrs > tt_os_dmet_yrs))



dft_dmet_surv %>% glimpse

# dft_dmet_surv %>%
#   filter(is.na(tt_os_dmet_yrs) | is.na(tt_cpt_dmet_yrs)) %>%
#   select(record_id, dx_to_dmets_yrs, tt_os_dmet_yrs, tt_os_dx_yrs)

# dft_dmet_surv %>% 
#   mutate(time_diff = tt_cpt_dmet_yrs - tt_os_dmet_yrs) %>%
#   ggplot(., aes(x = time_diff)) + stat_ecdf()


dft_dmet_surv %>% filter(tt_cpt_dmet_yrs > tt_os_dmet_yrs)
# THIS IS NOT RARE!

# Limitation of the method:
dft_dmet_surv %<>%
  filter(!(tt_cpt_dmet_yrs > tt_os_dmet_yrs))

# Create subsets based on the groups
dft_dmet_surv_trip_neg <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "Triple Negative")
dft_dmet_surv_hr_pos_her2_neg <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "HR+, HER2-")
dft_dmet_surv_her2_pos <- dft_dmet_surv %>%
  filter(bca_subtype_f_simple %in% "HER2+")


readr::write_rds(
  x = dft_dmet_surv,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet.rds')
)

readr::write_rds(
  x = dft_dmet_surv,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet.rds')
)







