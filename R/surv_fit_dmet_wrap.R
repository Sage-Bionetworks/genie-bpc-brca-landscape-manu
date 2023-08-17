
#' @description Big wrapper that creates the surv object, fits the survival,
#'   returns the datasets. Because the model is really specific to this one use
#'   case, we didn't put this wrapper into the /R directory.
#' @param dat A dataset with one row per {record_id}.
#' @param boot_rep The number of bootstrap resamples to do.
#' @param main_seed The seed used to draw the bootstrap seeds.
#' @param additional_features Any features to add to the X matrix.
surv_fit_dmet_wrap <- function(dat, boot_rep, main_seed, 
                               additional_features = character(0)) {
  y_dmet_os <- dat %>%
    select(
      tt_cpt_dmet_yrs_pos,
      tt_os_dmet_yrs,
      os_dx_status
    )
  
  x_dmet_os <- dat %>%
    # just to get the output in readable order we do an alphabet sort first:
    select(order(colnames(dat))) %>%
    select(
      matches("_mut$"),
      matches("_cna$"),
      matches("_fus$"),
      age_dx_c,
      stage_dx_iv_num,
      birth_year_c,
      white,
      hispanic,
      all_of(additional_features)
    ) %>%
    as.matrix(.)
  
  set.seed(main_seed) # for drawing each boot seed below
  dft_boot_dmet_os <- tibble(boot_ind = 1:boot_rep)%>%
    mutate(boot_seed = sample.int(n = 10^7, size = n()))
  
  rtn <- dft_boot_dmet_os %>%
    mutate(
      fit = purrr::map(
        .x = boot_seed,
        .f = (function(x) {
          fit_boot_cv_ltrc_lasso(
            x_mat = x_dmet_os,
            y_df = y_dmet_os,
            y_t = "tt_cpt_dmet_yrs_pos",
            y_t2 = "tt_os_dmet_yrs",
            y_e = "os_dx_status",
            boot_seed = x,
          )
        })
      )
    )
  
  return(rtn)
}