# Description: Creates folders for the analysis.

library(fs); library(here)

dir_create(here('data-raw', 'genomic'))
dir_create(here('data', 'genomic'))

fs::dir_create(here('data', 'survival', 'v2', 'prepared_data'))
fs::dir_create(here('data', 'survival', 'v2', 'fit_outputs', 'fit_summary'))

fs::dir_create(here('data', 'survival', 'age'))
fs::dir_create(here('data', 'survival', 'drug', 'simple_km_met'))

fs::dir_create(here('output', 'fig', 'manu'))
