library(magrittr)
library(dplyr)
library(janitor)

dft_cna <- readr::read_tsv(here('data-raw', 'data_CNA.txt'))
dft_cna %>% glimpse

dft_mut <- readr::read_tsv(here('data-raw', 'data_mutations_extended.txt'))
dft_mut %>% glimpse
dft_mut %>% tabyl(Hugo_Symbol) %>% arrange(desc(n)) %>% head(10)
# Probably start here.
dft_mut %>% glimpse

dft_fusion <- readr::read_tsv(here('data-raw', 'data_fusions.txt'))

dft_cpt <- readr::read_csv(
  here('data-raw', 'cancer_panel_test_level_dataset.csv')
)
