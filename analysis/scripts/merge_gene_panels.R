library(here)
library(purrr)
library(fs)
library(tibble)
library(dplyr)
library(yaml)

purrr::walk(fs::dir_ls('R'), .f = source)

vec_gene_panels <- fs::dir_ls('data-raw') %>%
  str_filter(., 'data_gene_panel_.*')

gp_all <- purrr::map_dfr(
  .x = vec_gene_panels,
  .f = tidy_gene_panel
)

# gp_all_wide <- gp_all %>%
#   tidyr::pivot_wider(
#     names_from = gene,
#     values_from = tested
#   )

saveRDS(
  object = gp_all,
  file = here('data', 'gene_panel_all.rds')
)
