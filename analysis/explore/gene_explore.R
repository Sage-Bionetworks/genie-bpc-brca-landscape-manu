library(magrittr)
library(dplyr)
library(janitor)
library(yaml)
library(here)
library(purrr)
library(fs)
library(ggplot2)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

dft_cpt <- readr::read_csv(here('data-raw', 'cancer_panel_test_level_dataset.csv'))
glimpse(dft_cpt)

dft_cna <- readr::read_tsv(here('data-raw', 'data_CNA.txt'))
dft_cna %>% glimpse

dft_mut <- readr::read_tsv(here('data-raw', 'data_mutations_extended.txt'))
dft_mut %>% glimpse
# some things which may help us filter a bit:
dft_mut$Consequence %>% unique
dft_mut$Variant_Classification %>% unique
dft_mut$Polyphen_Prediction %>% unique 
# Polyphen Score doesn't seem useful in addition due to...
ggplot(dft_mut, aes(x = Polyphen_Prediction, y = Polyphen_Score)) + geom_jitter()
dft_mut %>% tabyl(Hugo_Symbol) %>% arrange(desc(n)) %>% head(10)
# Probably start here.
dft_mut %>% glimpse

dft_fusion <- readr::read_tsv(here('data-raw', 'data_fusions.txt'))

dft_cpt <- readr::read_csv(
  here('data-raw', 'cancer_panel_test_level_dataset.csv')
)

dft_gene_mat <- readr::read_tsv(here('data-raw', 'data_gene_matrix.txt'))
dft_gene_mat %>% glimpse

dft_gp_test <- read_yaml(
  here('data-raw', 'data_gene_panel_DFCI-ONCOPANEL-1.txt')
)


tidy_gene_panel(
  here('data-raw', 'data_gene_panel_DFCI-ONCOPANEL-1.txt')
)





# Load some data to play with:
dft_pt <- readr::read_csv(here('data-raw', 'patient_level_dataset.csv'))
dft_ca_ind <- readr::read_csv(here('data-raw', 'cancer_level_dataset_index.csv'))
dft_dmet_timing <- get_dmet_timing(ca_ind_df = dft_ca_ind)
dft_gp_all <- readr::read_rds(here('data', 'gene_panel_all.rds'))


# Example:  Get all NGS which occur prior to metastatic cancer, or any which are
#   the first NGS test.
dft_cpt_dmet <- get_cpt_by_time(
  time_dat = dft_dmet_timing,
  time_var = "tt_y",
  cpt_dat = dft_cpt,
  always_keep_first = T
)

dft_cpt_dmet %>% glimpse
dft_cpt_dmet %>% tabyl(cpt_seq_assay_id)
