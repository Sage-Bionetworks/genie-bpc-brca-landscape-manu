# Cleans out files so we can start reproducibility from scratch

library(fs)

fs::dir_delete(here('data-raw'))

# selectively clear out the derived data - have some manual stuff in there still.
fs::dir_delete(here('data', 'clin_data_cohort'))
fs::dir_delete(here('data', 'genomic'))
fs::dir_delete(here('data', 'survival'))

fs::file_delete(here('data', 'gene_feat_long.rds'))
fs::file_delete(here('data', 'gene_feat_wide.rds'))
# not sure about this one:
fs::file_delete(here('data', 'drug_map.csv'))
