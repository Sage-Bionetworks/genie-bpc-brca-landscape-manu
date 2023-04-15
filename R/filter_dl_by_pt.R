#' @title Filter Data List by Progression Timing
filter_dl_by_pt <- function(d_list,
                            prog_timing) {
  # filter down to only the data you need to create plots of drug use.
  cohort_name <- names(data_list) 
  
  if (length(cohort_name) > 1) abort("Only designed for one cohort.")
  
  # filter
  data_list_exposed <- data_list[[cohort_name]]
  data_list_exposed <- data_list_exposed[c("pt_char",
                                           "ca_dx_index",
                                           "ca_drugs",
                                           "tumor_marker")]
  data_list_exposed[["pt_char"]] <- data_list_exposed[["pt_char"]] %>%
    filter(record_id %in% prog_timing$record_id)
  
  data_list_exposed[["ca_dx_index"]] <- data_list_exposed[["ca_dx_index"]] %>%
    filter(record_id %in% prog_timing$record_id)
  
  
  data_list_exposed[["ca_drugs"]] <- left_join(
    data_list_exposed[["ca_drugs"]],
    select(prog_timing, .data[["record_id"]], .data[["tt_d"]]),
    by = "record_id"
  ) %>%
    filter(record_id %in% prog_timing$record_id) %>%
    # only regimens that started after the time specified.
    filter(dx_reg_start_int >= tt_d) %>%
    select(-.data[["tt_d"]])
  
  
  data_list_exposed[["tumor_marker"]] <- left_join(
    data_list_exposed[["tumor_marker"]],
    select(prog_timing, .data[["record_id"]], .data[["tt_d"]]),
    by = "record_id"
  ) %>%
    filter(record_id %in% prog_timing$record_id) %>%
    # select only tumors with a date after progression:
    filter(dx_tm_days >= tt_d) %>%
    select(-.data[["tt_d"]])
  
  n_dat_list <- list()
  n_dat_list[[cohort_name]] <- data_list_exposed
  
  return(n_dat_list)
  
}

# filter_dl_by_pt(data_list,
#                 pt_pfs_i_or_m) %>%
#   lobstr::tree(., max_depth = 2)