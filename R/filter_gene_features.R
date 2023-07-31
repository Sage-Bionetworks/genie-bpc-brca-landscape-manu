#' @description Filters genes according to a sample list.  Further filters for nonzero variance after that, and a proportion of participants positive of at least 0.5% by default (arbitrary).
#' @param dat_gene_feat a "long" gene feature dataframe with columns cpt_genie_sample_id, feature, value.  Example features are "PTEN_mut or TP53_cna", values should be 0 or 1.
#' @param vec_samples a vector of samples
filter_gene_features <- function(dat_gene_feat, vec_samples, prop_filter = 0.005) {
  dat_gene_feat <- dat_gene_feat %>%
    filter(cpt_genie_sample_id %in% unique(vec_samples))
  
  # Require that there is some variance in the gene test results.
  dat_gene_feat %<>%   
    group_by(feature) %>%
    mutate(
      f_var = var(value, na.rm = T),
      prop_pos = mean(value %in% 1)
    ) %>%
    ungroup(.) %>%
    # 0.5% of dmet cases is about 3 people
    filter(f_var > 0 & prop_pos > 0.005)
  
  dat_gene_feat %<>% select(-c(f_var, prop_pos))
  
  return(dat_gene_feat)
}
