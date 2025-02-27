# Leveraging the work of my MSK colleagues to get the gene features.
# This takes the place of all the scripts that call oncogenic annotations, merge, etc.
library(purrr); library(here); library(fs);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_gp_all <- readr::read_rds(
  here('data', 'genomic', 'gene_panel_all.rds')
)
# hugo codes for genes covered in all panels
vec_hugo_ap <- dft_gp_all %>% 
  count(hugo) %>% 
  filter(n %in% max(n)) %>%
  pull(hugo)

# This is a base R style table with rownames.
msk_gam <- read.table(
  here('data', 'msk_manual_share', 'GAM_GENE_all_samples.txt')
)
msk_gam %<>%
  as_tibble(rownames = "sample_id")
msk_gam %<>%
  rename_with(
    .cols = everything(),
    .fn = \(colname) case_when(
       colname %in% 'CCND1.11q13.3' ~ 'CCND1',
       colname %in% 'NKX2.1' ~ 'NKX2-1',
       T ~ colname
    )
  )

msk_gam %<>%
  select(
    sample_id,
    contains("ERBB2"), # split into 3 - ok.
    all_of(vec_hugo_ap[!(vec_hugo_ap %in% "ERBB2")])
  )

# We expect 158 genes covered in all.
if (ncol(msk_gam) != 161) {
  cli_abort('Unexpected number of gene features (161 expected with ERBB2 split and sample_id included).  If the input data has changed then update or delete this check')
}

# For now I'll be renaming their features to match what I did previously.
# The reason for this is that I rely on feature names a bit in processing
#   the output and feature selection for models.
msk_gam %<>%
  rename_with(
    .cols = everything(),
    .fn = \(colname) {
      case_when(
        colname %in% "sample_id" ~ colname,
        colname %in% "ERBB2_MUT" ~ "ERBB2_mut",
        colname %in% "ERBB2_CNA" ~ "ERBB2_cna",
        colname %in% "ERBB2_FUS" ~ "ERBB2_fus",
        T ~ paste0(colname, "_gene")
      )
    }
  )

# There is, for some reason, a row which is 100% missing.  We can 
#   remove this by removing any feature, so let's do that and check.
nrow_before_filter <- nrow(msk_gam)
msk_gam %<>% filter(!is.na(ERBB2_mut)) 
nrow_after_filter <- nrow(msk_gam)

if (!anyNA(msk_gam)) {
  cli_alert_info("Gene matrix has no missing values after removing {nrow_before_filter - nrow_after_filter} rows.")
} else {
  cli_abort("Hacky row removal didn't work - fix.")
}

# There are still some slight differences from what I did before.
# I think I removed rows that were in removed people previously.
# Depending on the join used in subsequent scripts this might be fine.
        
readr::write_rds(
  x = msk_gam,
  file = here('data', 'genomic', 'gene_feat_oncogenic_msk.rds')
)

    