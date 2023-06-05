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
dft_cpt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_cpt.rds')
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

dft_mut <- left_join(
  dft_mut,
  select(dft_cpt, record_id, ca_seq, cpt_genie_sample_id, cpt_seq_assay_id),
  by = c(Tumor_Sample_Barcode = "cpt_genie_sample_id"),
  multiple = "error" # We expect one row per cpt_genie_sample_id.
)

dft_mut_long <- dft_mut %>%
  select(
    record_id, 
    ca_seq, 
    cpt_genie_sample_id = Tumor_Sample_Barcode, 
    hugo = Hugo_Symbol,
    cpt_seq_assay_id,
  )

# Just curious.  This should give the number of deleterious mutations per sample.
dft_mut_long %>% count(cpt_genie_sample_id, sort = T)


dft_mut %>% glimpse

# 
# dft_ca_ind <- readr::read_rds(
#   here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
# )
