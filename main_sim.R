# Description:  Top level workflow for the simulations.

library(purrr)
library(here)
library(fs)

purrr::walk(.x = fs::dir_ls('R'), .f = source)


# Simulation specific setup.  Does require some previous steps to be done.
source(here('analysis', 'script', 'sim_folder_setup.R'))
source(here('analysis', 'script', 'create_surv_data_real.R'))


# We declare a strategy for sessions here:
future::plan(strategy = multisession, workers = 6) 
# Cox univarate models (runs quickly - about 2 minutes with 6 workers)
source(here('analysis', 'script', 'run_method_univar_cox.R'))
source(here('analysis', 'script', 'eval_univar_cox.R'))


