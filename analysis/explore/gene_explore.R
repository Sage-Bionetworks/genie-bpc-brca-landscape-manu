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

# dft_gp_all <- readr::read_rds(
#   here('data', 'gene_panel_all.rds')
# )
# dft_gene_feat <- readr::read_rds(
#   here('data', 'gene_feat_long_imp.rds')
# )

dft_mut <- readr::read_rds(
  here('data', 'msk_box_derived', 'mut_reshape.rds')
)

dft_cna <- readr::read_rds(
  here('data', 'msk_box_derived', 'cna_reshape.rds')
)

rename_help <- function(x, suffix) {
  # need to do this so some of the stat functions don't freak out:
  x <- stringr::str_replace_all(x, "-", "_")
  x <- paste0(x, suffix)
  return(x)
}

# Give these names so we can tell the difference between CNA and mutation:
dft_mut %<>% 
  rename_at(
    .vars = vars(-stable_id),
    .funs = ~rename_help(., suffix = "_mut")
  ) 
     
dft_cna %<>% 
  rename_at(
    .vars = vars(-stable_id),
    .funs = ~rename_help(., suffix = "_cna")
  )    

dft_gene_feat <- full_join(dft_mut, dft_cna, by = "stable_id") %>%
  rename(cpt_genie_sample_id = stable_id)

# Go long to conform with previous code.
dft_gene_feat %<>%
  pivot_longer(
    cols = -cpt_genie_sample_id,
    names_to = "feature",
    values_to = "value"
  )
    

# dft_gene_feat %>% 
#   lapply(., FUN = (function(x) mean(is.na(x))))



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
  mutate(f_var = var(value, na.rm = T)) %>%
  ungroup(.) %>%
  filter(f_var > 0)

dft_gene_feat_dmet %<>% select(-f_var)

# Back to wide:
dft_gene_feat_dmet %<>%
  pivot_wider(
    names_from = "feature",
    values_from = "value",
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
count(dft_gene_feat_dmet, record_id, ca_seq, sort = T)
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
  full_join(., dft_bl, by = "record_id") 

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
  mutate(stage_dx_iv_num = if_else(stage_dx_iv %in% "Stage IV", 1, 0)) %>%
  select(
    PTEN_mut:WT1_cna,
    age_dx,
    stage_dx_iv,
    white,
    hispanic
  )





cv.glmnet(
  x = x_dmet,
  y = y_dmet,
  standardize = T,
  alpha = 1,
  family = "cox"
)

    