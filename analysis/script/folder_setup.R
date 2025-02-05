# Description: Creates folders for the analysis.

library(fs); library(here)

fs::dir_create(here("data-raw"), showWarnings = F)
dir.create(here('data', 'survival', 'prepared_data'), showWarnings = F)
dir.create(here('data', 'survival', 'fit_outputs'), showWarnings = F)
dir.create(here('data-raw', 'genomic'), showWarnings = F)
dir.create(here('data', 'genomic'), showWarnings = F)

fs::dir_create(here('data', 'survival', 'v2', 'prepared_data'))
fs::dir_create(here('data', 'survival', 'drug', 'fit_outputs'))
fs::dir_create(here('data', 'survival', 'drug', 'fit_outputs', 'fit_summary'))
fs::dir_create(here('data', 'survival', 'v2', 'fit_outputs', 'fit_summary'))
fs::dir_create(here('data', 'survival', 'age'))
fs::dir_create(here('data', 'survival', 'drug', 'simple_km_met'))

fs::dir_create(here('output', 'fig', 'manu'))
