# Description:  Setup simulation structure.

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

fs::dir_create(path = here("sim", "gen_data", "real_one"))
# set up results folders for each one:
fs::dir_create(path = here("sim", "results", "real_one", "cox_uni"))
fs::dir_create(path = here("sim", "results", "real_one", "reg_subsample"))
fs::dir_create(path = here("sim", "results", "real_one", "reg_boot"))

