library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

pt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_pt.rds')
)
ca_ind <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
)

tabyl(ca_ind, bca_subtype, ca_bca_her_summ)
# All people with borderline/equivocal/indeterminant are coded as NA (n=14)
# Same with people who are "Not applicable" and NA, test ordered by results not in chat, tests not done, unknown or no information.
tabyl(ca_ind, bca_subtype, ca_bca_er)
# Borderline is not included (n=1)
tabyl(ca_ind, bca_subtype, ca_bca_pr)
# HR positive = ER or PR, so borderlines are not included for the same person who was not ER positive.

ca_ind %>% 
  left_join(
    .,
    select(pt, record_id, birth_year, naaccr_race_code_primary, 
           naaccr_ethnicity_code, naaccr_sex_code),
    by = 'record_id'
  ) %>%
  select(
    birth_year, naaccr_race_code_primary, naaccr_ethnicity_code,
    naaccr_sex_code, age_dx, stage_dx, bca_subtype
  ) %>%
  mutate(bca_subtype = forcats::fct_na_value_to_level(bca_subtype, level = "Missing"),
         bca_st_miss = bca_subtype %in% "Missing") %>%
  tbl_summary(
    by = bca_st_miss
  )
