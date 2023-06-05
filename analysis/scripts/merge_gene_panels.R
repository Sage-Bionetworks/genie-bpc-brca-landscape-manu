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

saveRDS(
  object = gp_all,
  file = here('data', 'gene_panel_all.rds')
)




# Load in the mutation data and merge/process it.
cli::cli_alert_danger(
  "Mutation data processing was a guess to get the statistical steps set up.  
  Needs to be verified by someone with expertise!"
)
dft_cpt <- dft_cpt <- readr::read_csv(
  here('data-raw', 'cancer_panel_test_level_dataset.csv')
)
dft_mut <- readr::read_tsv(here('data-raw', 'data_mutations_extended.txt'))

dft_mut %<>%
  mutate(
    sift_flag = SIFT_Prediction %in% c("deleterious", "deleterious_low_confidence"),
    polyphen_flag = Polyphen_Prediction %in% c("probably_damaging", "possibly_damaging")
  ) %>%
  filter(sift_flag | polyphen_flag) %>%
  select(-c(sift_flag, polyphen_flag))

# The tumor sample barcodes seem to match up nicely with this column:
# dft_mut$Tumor_Sample_Barcode %in% dft_cpt$cpt_genie_sample_id %>% mean
dft_mut$Tumor_Sample_Barcode

left_join(
  dft_mut,
  select(dft_cpt, record_id, ca_seq, cpt_genie_sample_id),
  by = c(Tumor_Sample_Barcode = "cpt_genie_sample_id")
)

dft_cpt %>% filter(cpt_genie_sample_id %in% "GENIE-DFCI-093042-332548")
