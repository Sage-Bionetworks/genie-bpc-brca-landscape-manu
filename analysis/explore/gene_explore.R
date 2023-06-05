library(magrittr)
library(dplyr)
library(janitor)
library(yaml)
library(here)
library(purrr)
library(fs)
library(ggplot2)
library(tidyr)
library(survival)
library(glmnet)

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

dft_gp_all <- readr::read_rds(
  here('data', 'gene_panel_all.rds')
)
dft_gene_feat <- readr::read_rds(
  here('data', 'gene_feat_long.rds')
)




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

# The missingness and variance criteria build on McGough et. al. (2021), section
#. 4 titled "real-world data application"
# Further filter the genes down to those which were tested in at least 70% of panels.
dft_gene_feat_dmet %<>%
  group_by(hugo) %>%
  mutate(hugo_prop_tested = mean(tested %in% T)) %>%
  ungroup(.) %>%
  filter(hugo_prop_tested > 0.7)

cli::cli_alert_danger(
  text = "Filtering probably needs to be moved to post-merge."
)

# Also require that there is some variance in the gene test results.
dft_gene_feat_dmet %<>%   
  group_by(hugo) %>%
  mutate(hugo_var = var(variant, na.rm = T)) %>%
  ungroup(.) %>%
  filter(hugo_var > 0)

dft_gene_feat_dmet %<>% 
  select(
    -c(
      tested, 
      hugo_prop_tested, 
      hugo_var,
      cpt_seq_assay_id
    )
  ) 

dft_gene_feat_dmet %<>%
  pivot_wider(
    names_from = "hugo",
    values_from = "variant",
    values_fill = NA
  )
# Should be all 1's:
# count(dft_gene_feat_dmet, cpt_genie_sample_id, sort = T)
# count(dft_cpt_dmet, cpt_genie_sample_id, sort = T)

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
  select(record_id, white, hispanic)
      

dft_bl <- dft_ca_ind %>%
  select(
    record_id, 
    ca_seq, # works because we've already selected one row per person.
    age_dx, 
    stage_dx_iv, 
    os_dx_status, 
    tt_os_dx_yrs
  ) %>%
  full_join(., dft_dmet_bl, by = "record_id") 

dft_dmet_surv <- left_join(
  dft_gene_feat_dmet,
  dft_bl,
  by = c("record_id", "ca_seq"),
  multiple = "error"
)

# Limitation of the method:
dft_dmet_surv %<>%
  filter(!(dx_cpt_rep_yrs > tt_os_dx_yrs))






y_dmet <- with(
  dft_dmet_surv,
  Surv(
    time = dx_cpt_rep_yrs,
    time2 = tt_os_dx_yrs,
    event = os_dx_status
  )
)

x_dmet <- dft_dmet_surv %>% 
  select(
    ABL1:WT1,
    age_dx,
    stage_dx_iv,
    white,
    hispanic
  )


cli::cli_alert_danger(text = "Need to do imputation still!!!!")



cv.glmnet(
  x = x_dmet,
  y = y_dmet,
  standardize = T,
  alpha = 1,
  family = "cox"
)

    