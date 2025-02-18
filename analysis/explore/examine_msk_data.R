
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

dir_input <- here('data', 'msk_manual_share')
dir_out <- here('/Users/apaynter/main/projects/genie/2024_12_11_breast_manu/compare_code')

gam_gene <- read.table(
  here(dir_input, 'GAM_GENE_all_samples.txt')
)

# probably don't need this one but here it is:
# gam_pathway <- read.table(
#   here(dir_input, 'GAM_PATHWAY_all_samples.txt')
# )

clin_samp <- readr::read_tsv(
  here(dir_input, 'master_clinical_table_all_samples.tsv')
)

included_in_landscape <- clin_samp %>%
  # always remember to code T/F columns as "Keep" and "Exclude"...
  filter(Keep_for_Landscape %in% "Keep") %>%
  select(PATIENT_ID, SAMPLE_ID)

if (included_in_landscape$PATIENT_ID %>% duplicated %>% any) {
  cli_abort("Duplicate PATIENT_IDs in included_in_landscape")
}


# do a bit of processing on the genomics - finding samples with zero alts.
gam_gene %<>%
  as_tibble(rownames = "sample_id")
gam_gene_long <- gam_gene %>% 
  pivot_longer(
    cols = -sample_id,
    names_to = "feature",
    values_to = "value"
  )
samp_gene_ct <- gam_gene_long %>%
  group_by(sample_id) %>%
  summarize(n_alt = sum(value, na.rm = T))
flat_cases <- samp_gene_ct %>%
  filter(n_alt %in% 0) %>%
  pull(sample_id)



clin_samp <- clin_samp %>%
  mutate(
    pt_lvl_keep = PATIENT_ID %in% included_in_landscape,
    is_flat_case = SAMPLE_ID %in% flat_cases,
    bca_subtype_missing = is.na(BCA_SUBTYPE),
    is_sarcoma = CANCER_TYPE %in% "Breast Sarcoma"
  ) %>%
  left_join(., samp_gene_ct, by = c(SAMPLE_ID = "sample_id"))

excluded_unexplained <- clin_samp %>%
  filter(
    !pt_lvl_keep & 
      !is_flat_case &
      !bca_subtype_missing &
      !is_sarcoma & 
      Keep_for_Landscape %in% "Exclude"
  )


excluded_unexplained %>%
  rename_all(tolower) %>%
  select(patient_id, sample_id, cancer_type, sex, bca_subtype, n_alt) %>%
  View(.)

  # filter(Keep_for_Landscape %in% "Exclude" & 
  #          !(PATIENT_ID %in% included_in_landscape)) %>%
  # select(PATIENT_ID, SAMPLE_ID)

included_in_landscape %<>% rename_all(tolower)
excluded_with_no_duplicate %<>% rename_all(tolower)


readr::write_rds(
  included_in_landscape,
  here(dir_out, 'msk_included_in_landscape.rds')
)
readr::write_rds(
  excluded_unexplained,
  here(dir_out, 'msk_excluded_pts.rds')
)




clin_samp %>%
  mutate(pt_has_another_sample = PATIENT_ID %in% included_in_landscape)
  filter(Keep_for_Landscape %in% "Exclude" & 
           !(PATIENT_ID %in% pts_in_landscape)) %>%
  select(PATIENT_ID, SAMPLE_ID)


gam_gene %>%
  summarize(
    across(
      .cols = everything(),
      .fns = \(z) mean(is.na(z))
    )
  ) 



gp_all <- readr::read_rds(
  here('data', 'genomic', 'alt_test_full.rds')
)
alt_test <- readr::read_rds(
  here('data', 'genomic', 'alt_test_full.rds')
)
# I'm not sure why it has duplicates, but we can fix that:
alt_test %<>%
  filter(alt_type %in% "Mutation") %>%
  group_by(sample_id, hugo) %>%
  slice(1) %>%
  ungroup(.)

gam_gene_long %>%
  filter(sample_id %in% included_in_landscape$sample_id) %>%
  group_by(feature) %>%
  summarize(
    n_samp_tested = sum(!is.na(value))
  ) %>%
  filter(n_samp_tested == length(unique(included_in_landscape$sample_id))) %>% View(.)

test_compare <- gam_gene_long %>%
  filter(sample_id %in% included_in_landscape$sample_id) %>%
  mutate(
    tested_msk = if_else(is.na(value), F, T)
  ) %>%
  full_join(
    ., 
    (alt_test %>%
       filter(sample_id %in% included_in_landscape$sample_id) %>%
       select(sample_id, feature = hugo, tested_alex = tested,
              cpt_seq_assay_id)),
    by = c("sample_id", "feature")
  )

test_compare %>%
  group_by(feature) %>%
  mutate(across(.cols = c(tested_msk, tested_alex),
                .fns = \(x) if_else(is.na(x), F, x))) %>%
  summarize(
    pct_test_msk = mean(tested_msk),
    pct_agree = mean(tested_msk == tested_alex)
  ) %>%
  arrange(pct_agree) %>%
  View(.)

