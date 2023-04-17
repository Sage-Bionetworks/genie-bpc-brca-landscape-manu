# Description: A sandbox for looking at and playing with the datasets.
#   This script will never be a part of the final analysis workflow.
# Author: Alex Paynter

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
library(lobstr) # just for lobstr::tree - more readable playing

purrr::walk(.x = here("R", dir(here("R"))), .f = source)

data_list <- readr::read_rds(here("data-raw", "data_list.rds"))

dft_pt_char <- data_list$BrCa_v1.2$pt_char
dft_ca_dx_index <- data_list$BrCa_v1.2$ca_dx_index
dft_ca_dx_non_index <- data_list$BrCa_v1.2$ca_dx_non_index
dft_ca_drugs <- data_list$BrCa_v1.2$ca_drugs

dft_prissmm_imaging <- data_list$BrCa_v1.2$prissmm_imaging
dft_prissmm_pathology <- data_list$BrCa_v1.2$prissmm_pathology
dft_prissmm_md <- data_list$BrCa_v1.2$prissmm_md
dft_tumor_marker <- data_list$BrCa_v1.2$tumor_marker








# Try something as a fix for sunburst errors:
nrow(data_list$BrCa_v1.2$ca_dx_index)

data_list$BrCa_v1.2$ca_dx_index <- data_list$BrCa_v1.2$ca_dx_index %>%
  group_by(record_id) %>% 
  slice(1) %>%
  ungroup




pt_pfs_i_or_m <- get_progressed_timing(ca_dat = dft_ca_dx_index, 
                                       prefix = "pfs_m_adv")
faux_cohort <- filter_dl_by_pt(data_list,
                               pt_pfs_i_or_m) 


# Closer, at least it runs:
test <- drug_regimen_sunburst(
  data_list$BrCa_v1.2[c("pt_char", "ca_dx_index", "ca_drugs", "cpt")],
  faux_cohort,
  max_n_regimens = 2
)

# The real, valuable output of the function is this:
test$treatment_history

# The sunplot that they return can be reproduced in one line from there,
test$sunburst_plot
test$treatment_history %>% sunburstR::sunburst(., legend = T)
# which is great, because we can probably jazz it up from here a bit:
set.seed(12) # for color sampling.
sunburstR::sunburst(
  test$treatment_history,
  legend = T,
  colors = sample(
    viridisLite::turbo(n = nrow(test$treatment_history),
                         begin = 0, end = 1))
)

