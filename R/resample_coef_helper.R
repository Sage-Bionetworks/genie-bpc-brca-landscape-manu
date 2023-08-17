resample_coef_helper <- function(
    dat, 
    col_with_fits = "fit",
    filter_log_hr_above = 4.60517, # HR = 100 
    stab_thresh = NULL,
    exp_coefs = T
) {
  dat %<>%
    mutate(
      coef_mat = purrr::map(
        .x = .data[[col_with_fits]], 
        .f = get_cv_lasso_coefs)
    ) %>%
    select(-fit) %>%
    unnest(coef_mat)
  
  dat %<>%
    mutate(
      feature = forcats::fct_inorder(feature),
    ) %>%
    filter(log_hr < filter_log_hr_above & log_hr > -filter_log_hr_above) %>%
    group_by(feature) %>%
    summarize(
      stability = mean(abs(log_hr) > 0.01),
      log_hr = mean(log_hr),
      n = n()
    ) %>%
    mutate(
      hr = exp(log_hr)
    )
  
  if (exp_coefs) {
    dat %<>% select(-log_hr)
  } else {
    dat %<>% select(-hr)
  }
  
  if (!is.null(stab_thresh)) {
    dat %<>% 
      filter(stability >= stab_thresh)
  }
  
  return(dat)
}
