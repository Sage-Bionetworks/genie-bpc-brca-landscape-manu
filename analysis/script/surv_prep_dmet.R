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

# Doing this now for my cohort tracking - possibly could do it at the start.
dft_ca_ind %<>% 
  mutate(
    bca_subtype_f_simple = forcats::fct_na_value_to_level(
      f = bca_subtype_f_simple,
      level = "(NC or NR)"
    )
  )
    
    
surv_cohort_track_help <- function(dat, state, var = "bca_subtype_f_simple") {
  dat %>%
    tabyl(all_of(var)) %>%
    adorn_totals(., name = "All") %>%
    select(1:2) %>%
    mutate(
      state = state
    )
}

dft_surv_consort <- surv_cohort_track_help(dft_ca_ind, state = "start")


# dft_cpt %>%
#   select(record_id, cpt_order_int, cpt_report_int) %>%
#   filter(is.na(cpt_order_int) | is.na(cpt_report_int))
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

dft_clin_char %<>% 
  mutate(
    dx_to_dmets_yrs = if_else(stage_dx_iv %in% "Stage IV", 0, dx_to_dmets_yrs),
    tt_os_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ tt_os_dx_yrs - dx_to_dmets_yrs,
    )
    # pfs is already relative to stage IV or dmet date.
  ) 

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




# Do the data processing for all people first:
dft_gene_comb_all <- combine_cpt_gene_feat(
  dat_gene_feat = dft_gene_feat_dmet,
  dat_cpt = dft_cpt_dmet
)

dft_dmet_surv_all <- combine_clin_gene(
  dat_gene = dft_gene_comb_all,
  dat_clin = dft_clin_char
) %>%
  add_specific_dmet_vars(.)

dft_surv_consort <- bind_rows(
  dft_surv_consort,
  surv_cohort_track_help(dat = dft_dmet_surv_all, state = "dmet")
)

dft_dmet_surv_all %<>%
  filter(!(tt_cpt_dmet_yrs > tt_os_dmet_yrs))

dft_surv_consort <- bind_rows(
  dft_surv_consort,
  surv_cohort_track_help(dat = dft_dmet_surv_all, state = "dmet, CPT <= OS follow-up")
) 

dft_surv_consort %<>%
  mutate(state = forcats::fct_inorder(f = state)) %>%
  select(bca_subtype_f_simple, state, n) %>%
  arrange(bca_subtype_f_simple, state)

readr::write_rds(
  file = here('data', 'survival', 'prepared_data', 'surv_consort_dmet.rds'),
  x = dft_surv_consort
)
  

#' @description Filters gene features for non-zero variance (obvious criterion)
#'   and so that a certain proportion of participants (0.5% by default,
#'   arbitrary) are positive.
#' @param dat a dataframe with one row per case ({record_id, ca_seq})
#' @param prop_filter the proportion of participants which must be positive for
#'   the feature to stay in the dataset.
#' @param additional_features By default, the function will apply this filtering
#'   to any column which ends in '_mut', '_cna' or '_fus', a convention for the
#'   types of genomic alterations we're considering.  If there are additional
#'   custom gene features to be considered, add their names here.
filter_gene_features <- function(dat, prop_filter = 0.005, additional_features = character(0)) {

  
  dat_gene <- dat %>%
    select(
      record_id, ca_seq,
      matches("_mut$"),
      matches("_cna$"),
      matches("_fus$"),
      all_of(additional_features),
    ) 
  
  # pull out the non-gene features, we can combine them again at the end.
  cols_not_gene <- setdiff(names(dat), names(dat_gene)) 
  dat_not_gene <- select(dat, record_id, ca_seq, all_of(cols_not_gene))
  
  
  dat_gene %<>%
    pivot_longer(
      cols = -c(record_id, ca_seq),
      names_to = "feature",
      values_to = "value"
    )
  
  # Require that there is some variance in the gene test results.
  dat_gene %<>%   
    group_by(feature) %>%
    mutate(
      f_var = var(value, na.rm = T),
      prop_pos = mean(value %in% 1)
    ) %>%
    ungroup(.) %>%
    # 0.5% of dmet cases is about 3 people
    filter(f_var > 0 & prop_pos > prop_filter)
  
  dat_gene %<>% select(-c(f_var, prop_pos))
  
  dat_gene <- pivot_wider(
    dat_gene,
    names_from = "feature", values_from = "value"
  )
  
  rtn <- left_join(
    dat_not_gene,
    dat_gene,
    by = c("record_id", "ca_seq")
  )
  
  return(rtn)
}



# Before we filter the genes for variance/proportions, make the subgroup datasets:
dft_dmet_surv_trip_neg <- dft_dmet_surv_all %>%
  filter(bca_subtype_f_simple %in% "Triple Negative") %>%
  filter_gene_features(.)

dft_dmet_surv_hr_pos_her2_neg <- dft_dmet_surv_all %>%
  filter(bca_subtype_f_simple %in% "HR+, HER2-") %>%
  filter_gene_features(.)

dft_dmet_surv_her2_pos <- dft_dmet_surv_all %>%
  filter(bca_subtype_f_simple %in% "HER2+") %>%
  filter_gene_features(.)

# Now we can filter the overall data as well:
dft_dmet_surv_all %<>% filter_gene_features(.)

# Note: If you decided to require both 0.5% in the overall cohort and 0.5% in the
#   subgroup, you could move the filter on "_all" to be before the subgroup ones.

 

# Save the datasets.
readr::write_rds(
  x = dft_dmet_surv_all,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet_all.rds')
)
readr::write_rds(
  x = dft_dmet_surv_trip_neg,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet_trip_neg.rds')
)
readr::write_rds(
  x = dft_dmet_surv_hr_pos_her2_neg,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet_hr_pos_her2_neg.rds')
)
readr::write_rds(
  x = dft_dmet_surv_her2_pos,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet_her2_pos.rds')
)

# Comments on revised process that applies the genetic filters to each dataset separately:
# Before - after numbers on each dataset for columns:
# - all: 107 - 104
# - trip_neg: 107 - 153
# - hr_pos_her2_neg: 107 - 95
# - her2_pos: 107 - 152

# As it should be, the row count was the same before and after, we just lose/gain some features.







