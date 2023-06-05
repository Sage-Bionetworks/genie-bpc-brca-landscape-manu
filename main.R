# Description:  Top level workflow for the project.  Eventually this should
#  be automated, clean.  For now it's a list of comments and source files.

library(purrr)
library(here)
library(fs)

# Load all helper functions
purrr::walk(.x = fs::dir_ls('R'), .f = source)

source(here('analysis', 'scripts', 'get_raw_data.R'))
source(here('analysis', 'scripts', 'filter_data_for_cohort.R'))
