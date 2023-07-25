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

dft_cpt %>%
  select(record_id, cpt_order_int, cpt_report_int) %>%
  mutate(day_dist = cpt_order_int - cpt_report_int) %>%
  pull(day_dist) %>%
  summary






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

dft_gene_feat_dmet <- dft_gene_feat %>%
  filter(cpt_genie_sample_id %in% unique(dft_cpt_dmet$cpt_genie_sample_id))

# Require that there is some variance in the gene test results.
dft_gene_feat_dmet %<>%   
  group_by(feature) %>%
  mutate(
    f_var = var(value, na.rm = T),
    prop_pos = mean(value %in% 1)
  ) %>%
  ungroup(.) %>%
  # 0.5% of dmet cases is about 3 people
  filter(f_var > 0 & prop_pos > 0.005)

dft_gene_feat_dmet %<>% select(-c(f_var, prop_pos))

# Take this chance to arrange them logically:
# dft_gene_feat_dmet %>%
#   arrange(feature) %>%
#   arrange(desc(str_detect(feature, "_mut$") * 3 + 
#                  str_detect(feature, "_cna$") * 2 + 
#                  str_detect(feature, "_fus$") * 1)
#   )
# Did not work for some reason.

# Back to wide:
dft_gene_feat_dmet %<>%
  pivot_wider(
    names_from = "feature",
    values_from = "value",
    values_fill = NA
  )
# Should be all 1's:
if (max(pull(count(dft_gene_feat_dmet, cpt_genie_sample_id, sort = T), n)) > 1) {
  cli::cli_abort("Duplicate gene features after wide pivot")
}
if (max(pull(count(dft_cpt_dmet, cpt_genie_sample_id, sort = T), n)) > 1) {
  cli::cli_abort("Duplicate NGS rows.")
}
  

dft_gene_feat_dmet <- left_join(
  dft_cpt_dmet,
  dft_gene_feat_dmet,
  by = c("cpt_genie_sample_id")
)

# Now we have to deal with people that have more than one test:
# count(dft_gene_feat_dmet, record_id, ca_seq, sort = T)
# OK, we have two people, so it doesn't seem too important here.  
# Let's drop their second tests.
dft_gene_feat_dmet %<>%
  filter(is_first_cpt) 
# Problem solved - verify with most recent count command.



# Could add other covariates here from dft_pt, but just to start:
dft_bl <- dft_pt %>%
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


dft_bl <- dft_ca_ind %>%
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
  full_join(., dft_bl, by = "record_id") 

dft_bl %>% mutate(na_check = is.na(dx_to_dmets_yrs)) %>% tabyl(na_check)

dft_bl %<>% 
  mutate(
    dx_to_dmets_yrs = if_else(stage_dx_iv %in% "Stage IV", 0, dx_to_dmets_yrs),
    tt_os_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ tt_os_dx_yrs - dx_to_dmets_yrs,
    )
    # pfs is already relative to stage IV or dmet date.
  ) 
dft_dmet_surv <- left_join(
  dft_gene_feat_dmet,
  dft_bl,
  by = c("record_id", "ca_seq"),
  multiple = "error"
)

# cleaning up a few of these for later interpretation:
dft_dmet_surv <- dft_dmet_surv %>%
  mutate(
    tt_cpt_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ dx_cpt_rep_yrs - dx_to_dmets_yrs,
    ),
    # glmnet won't take negative times:
    tt_cpt_dmet_yrs_pos = if_else(tt_cpt_dmet_yrs < 0, 0, tt_cpt_dmet_yrs),
    stage_dx_iv_num = if_else(stage_dx_iv %in% "Stage IV", 1, 0),
    age_dx_c = age_dx - 40, # approximately centered.
    birth_year_c = birth_year - 1970, # approximately centered.
  )

dft_dmet_surv %>%
  filter(is.na(tt_os_dmet_yrs) | is.na(tt_cpt_dmet_yrs)) %>%
  select(record_id, dx_to_dmets_yrs, tt_os_dmet_yrs, tt_os_dx_yrs)

dft_dmet_surv %>% 
  mutate(time_diff = tt_cpt_dmet_yrs - tt_os_dmet_yrs) %>%
  ggplot(., aes(x = time_diff)) + stat_ecdf()

# Limitation of the method:
dft_dmet_surv %<>%
  filter(!(tt_cpt_dmet_yrs > tt_os_dmet_yrs))


readr::write_rds(
  x = dft_dmet_surv,
  file = here('data', 'survival', 'prepared_data', 'surv_dmet.rds')
)









# dft_dmet_timing %>% 
#   mutate(developed_dmet = 1) %>% 
#   select(record_id, developed_dmet, dmet_yrs) %>%
#   left_join(
#     (select(dft_ca_ind, record_id, ca_seq, bca_subtype_f_simple, 
#             tt_os_dmet_yrs) %>%
#        tt_os_dmet_yrs = case_when(
#          is.na(dx_to_dmets_yrs) ~ NA_real_,
#          T ~ tt_os_dx_yrs - dx_to_dmets_yrs,
#        )
#     ),
#     .,
#     by = "record_id"
#   )
# # Getting the numbers for a consort-like diagram here.
# vec_surv_consort_dmet_all <- character(0L)
# vec_surv_consort_dmet_all["Cohort"] <-  dft_ca_ind
