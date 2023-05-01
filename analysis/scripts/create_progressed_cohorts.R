library(cli) # prevents an error on version reliance.
library(readr)
library(vctrs)
library(rlang)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(janitor)
library(glue)
library(genieBPC)
library(ggplot2)
library(sunburstR)
library(huxtable)

purrr::walk(.x = here("R", dir(here("R"))), .f = source)


data_list <- readr::read_rds(here("data-raw", "data_list.rds"))

dft_med_classified <- readr::read_csv(
  file = here("data", "tx_by_class_metastatic_filt_cohort.csv")
)

dft_med_classified <- dft_med_classified %>%
  select(
    record_id = PATIENT_ID,
    regimen_number = REGIMEN_NUMBER,
    ca_seq = CA_SEQ,
    agent = AGENT,
    regimen = REGIMEN,
    class1 = Class1,
    class1.1 = Class1.1,
    exclude = Exclude
  )

dft_med_classified %>% pull(regimen) %>% unique %>% length

dft_reg_map <- dft_med_classified %>%
  group_by(record_id, regimen_number, ca_seq, regimen) %>%
  summarize(
    regimen_c1 = paste0(sort(unique(class1)), collapse = ", "),
    regimen_c1.1 = paste0(sort(unique(class1.1)), collapse = ", "),
    .groups = "drop"
  ) %>%
  group_by(regimen) %>% 
  slice(1) %>%
  select(regimen, regimen_c1, regimen_c1.1) %>%
  ungroup()

# The sunburst code uses dashes as special characters so "pd-l1" and "ctla-4" will cause problems.  We can just remove those without too much loss of info:
dft_reg_map %<>%
  mutate(regimen_c1.1 = stringr::str_replace(regimen_c1.1,
                                             "-",
                                             "")) 

dft_drug_map <- dft_med_classified %>%
  select(agent, class1, class1.1) %>%
  mutate(class1.1 = stringr::str_replace(class1.1,
                                             "-",
                                             "")) %>%
  group_by(agent) %>% 
  slice(1) %>%
  ungroup()



data_list$BrCa_v1.2$ca_dx_index <- data_list$BrCa_v1.2$ca_dx_index %>%
  group_by(record_id) %>% 
  slice(1) %>%
  ungroup

data_list_c1 <- data_list
data_list_c1.1 <- data_list

# replace the list of drugs with the list of regimens from this set:
data_list_c1$BrCa_v1.2$ca_drugs <- data_list_c1$BrCa_v1.2$ca_drugs %>%
  mutate(regimen_drugs = as.character(regimen_drugs)) %>%
  left_join(
    .,
    select(dft_reg_map, regimen, regimen_c1),
    by = c(regimen_drugs = "regimen")
  ) %>%
  # Can't filter here due to the genieBPC function below breaking:
  # mutate(regimen_c1 = if_else(is.na(regimen_c1) | regimen_c1 %in% "",
  #                             "leuprolide, inv_drug, other",
  #                             regimen_c1)) %>%
  filter(!(regimen_c1 %in% "") & !is.na(regimen_c1)) %>%
  mutate(regimen_drugs = factor(regimen_c1)) %>%
  select(-regimen_c1)



# Fix the regimen numbering - required for the genieBPC package function:
data_list_c1$BrCa_v1.2$ca_drugs <- data_list_c1$BrCa_v1.2$ca_drugs %>%
  arrange(record_id, ca_seq, regimen_number) %>%
  group_by(record_id, ca_seq) %>%
  mutate(regimen_number = 1:n()) %>%
  ungroup


data_list_c1$BrCa_v1.2$ca_drugs <- 
  data_list_c1$BrCa_v1.2$ca_drugs %>% 
  mutate(
    across(
      .cols = matches("^drugs_drug_[0-9]$"),
      .fns = (function(x) {
        x <- as.character(x)
        x <- stringr::str_replace(x, "\\(.*", "")
        x <- map_drug(x, dft_drug_map, type = "class1")
        return(x)
      })
    )) 







# Repeat for c1.1:
data_list_c1.1$BrCa_v1.2$ca_drugs <- data_list_c1.1$BrCa_v1.2$ca_drugs %>%
  mutate(regimen_drugs = as.character(regimen_drugs)) %>%
  left_join(
    .,
    select(dft_reg_map, regimen, regimen_c1.1),
    by = c(regimen_drugs = "regimen")
  ) %>%
  # It does also work to do this, or something similar, if you want to include Leuprolide or Cyclophosphamide later on.
  # mutate(regimen_c1 = if_else(is.na(regimen_c1) | regimen_c1 %in% "",
  #                             "inv drug or leup", 
  #                             regimen_c1)) %>%
  filter(!(regimen_c1.1 %in% "") & !is.na(regimen_c1.1)) %>%
  mutate(regimen_drugs = factor(regimen_c1.1)) %>%
  select(-regimen_c1.1)

data_list_c1.1$BrCa_v1.2$ca_drugs <- data_list_c1.1$BrCa_v1.2$ca_drugs %>%
  arrange(record_id, ca_seq, regimen_number) %>%
  group_by(record_id, ca_seq) %>%
  mutate(regimen_number = 1:n()) %>%
  ungroup


data_list_c1.1$BrCa_v1.2$ca_drugs <- 
  data_list_c1.1$BrCa_v1.2$ca_drugs %>% 
  mutate(
    across(
      .cols = matches("^drugs_drug_[0-9]$"),
      .fns = (function(x) {
        x <- as.character(x)
        x <- stringr::str_replace(x, "\\(.*", "")
        x <- map_drug(x, dft_drug_map, type = "class1.1")
        return(x)
      })
    )) 





dft_dmet_timing <- get_dmet_timing(
  ca_ind_df = data_list_c1$BrCa_v1.2$ca_dx_index
)
prog_cohort_c1 <- filter_dl_by_pt(
  d_list = data_list_c1,
  dft_dmet_timing
)

dft_dmet_timing_c1.1 <- get_dmet_timing(
  ca_ind_df = data_list_c1.1$BrCa_v1.2$ca_dx_index
)


data_list_c1.1$BrCa_v1.2$ca_drug
prog_cohort_c1.1 <- filter_dl_by_pt(
  d_list = data_list_c1.1,
  dft_dmet_timing_c1.1
)

save(
  data_list,
  data_list_c1,
  data_list_c1.1,
  dft_dmet_timing,
  dft_dmet_timing_c1.1,
  dft_med_classified,
  dft_drug_map,
  dft_reg_map,
  prog_cohort_c1,
  prog_cohort_c1.1,
  file = here("data", "prog_cohorts.Rda")
)



