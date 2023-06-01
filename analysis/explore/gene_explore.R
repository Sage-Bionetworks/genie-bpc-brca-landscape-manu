library(magrittr)
library(dplyr)
library(janitor)

dft_

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

dft_gene_mat <- readr::read_tsv(here('data-raw', 'data_gene_matrix.txt'))
dft_gene_mat %>% glimpse

dft_gp_test <- readLines(
  here('data-raw', 'data_gene_panel_DFCI-ONCOPANEL-1.txt')
) 
