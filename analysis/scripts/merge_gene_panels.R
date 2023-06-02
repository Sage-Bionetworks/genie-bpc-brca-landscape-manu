library(here)
library(purrr)
library(fs)
library(tibble)
library(dplyr)

purrr::walk(fs::dir_ls('R'), .f = source)

dft_gene_panels <- fs::dir_ls('data-raw') %>%
  str_filter(., 'data_gene_panel_.*') %>%
  tibble(path = .)

dft_gene_panels %>% slice(1) %>% pull(path) %>% tidy_gene_panel(.)

