library(purrr); library(here); library(fs);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)


# dft_ca_ind <- readr::read_rds(
#   here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
# )
dft_cpt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_cpt.rds')
)
dft_gene_feat_wide <- readr::read_rds(
  here('data', 'genomic', 'gene_feat_oncogenic.rds')
)
# Clinical characteristics already processed some:
dft_clin_char <- readr::read_rds(
  here(
    'data', 'survival', 'v2', 'prepared_data',
    'clin_char.rds'
  )
)

dft_drug_feas_surv <- readr::read_rds(
  here('data', 'clin_data_cohort', 'drug_feas_surv.rds')
)



# Do we actually need this dataset?
# dft_ca_ind %<>%
#   mutate(
#     bca_subtype_f_simple = forcats::fct_na_value_to_level(
#       f = bca_subtype_f_simple,
#       level = "(NC or NR)"
#     )
#   )

dft_clin_char %<>%
  # take out some variables that could lead to confusion:
  select(
    -contains('os_'), 
    -contains("pfs_"),
    -bca_subtype, -bca_subtype_f,
    -dx_to_dmets_yrs, # already have this in drug feas dataset.
    -ca_seq # technially not a problem, but every record has just one ID now.
  ) %>%
  mutate(
    bca_subtype_f_simple = forcats::fct_na_value_to_level(
      f = bca_subtype_f_simple,
      level = "(NC or NR)"
    )
  )



dft_drug_surv <- dft_drug_feas_surv %>% 
  filter(crit_all_os) %>%
  select(-contains("crit"))

dft_drug_surv <- left_join(
  dft_drug_surv,
  dft_clin_char,
  by = "record_id"
)

dft_drug_surv %>%
  mutate(
    # this age is not exact.  It's a best guess since age_dx is rounded.
    age_drug_start = age_dx + drug_dx_start_int_yrs
  ) %>%
  # just to avoid mistakes:
  select(-age_dx)


  











# POC: do the CDK inhibitor case.
dft_drug_feas_cdk <- dft_drug_feas_surv %>%
  filter(class_comp %in% 'CDK inhibitor') %>%
  # participants who meet all time-based criteria for the OS analysis:
  filter(crit_all_os) %>%
  select(record_id, drug_dx_start_int_yrs)
# Function below requires ca_seq, which we can add back in:
dft_drug_feas_cdk %<>%
  left_join(
    .,
    select(dft_ca_ind, record_id, ca_seq),
    by = c("record_id")
  )
# Get all NGS which occur prior to drug use.
dft_cpt_drug <- get_cpt_by_time(
  time_dat = dft_drug_feas_cdk,
  time_var = "drug_dx_start_int_yrs",
  cpt_dat = dft_cpt,
  always_keep_first = F
)

dft_drug_feas_cdk %<>% select(-ca_seq)

dft_cpt_drug %<>%
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
    # skipping for now:
    # pfs_i_and_m_adv_status,
    # tt_pfs_i_and_m_adv_yrs
  ) %>%
  full_join(., dft_clin_char, by = "record_id") 

# drug specific:
dft_clin_char <-
  left_join(
    dft_drug_feas_cdk,
    dft_clin_char,
    by = "record_id"
  )

dft_clin_char %<>% 
  mutate(
    tt_os_drug_yrs = tt_os_dx_yrs - drug_dx_start_int_yrs
  ) %>%
  select(-tt_os_dx_yrs) # avoid mistakes
    




# For the analyses of the drug it makes the most sense to me to keep the most recent test:
dft_gene_comb_drug <- combine_cpt_gene_feat(
  dat_gene_feat = dft_gene_feat,
  dat_cpt = dft_cpt_drug,
  keep_only_first = F
) %>%
  group_by(record_id) %>%
  arrange(dx_cpt_rep_yrs) %>%
  slice(n()) %>%
  ungroup

dft_dmet_surv_drug <- combine_clin_gene(
  dat_gene = dft_gene_comb_drug,
  dat_clin = dft_clin_char
) 

# We've already done all the time filtering at this point, no need to do it again.




# Before we filter the genes for variance/proportions, make the subgroup datasets:
dft_dmet_surv_drug %<>% 
  filter_gene_features(., 0.02)


readr::write_rds(
  x = dft_dmet_surv_drug,
  file = here('data', 'survival', 'drug', 'prepared_data', 'dmet_surv_cdk.rds')
)
  








