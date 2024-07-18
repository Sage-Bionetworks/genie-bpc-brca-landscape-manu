# Description:  This was a request from the breast cancer group as a compliment
#  to Evan's analysis.  I did not build on surv_prep_dmet because we don't 
#  need any of the complexity of genomic data here.

library(here); library(purrr); library(fs)
purrr::walk(.x = fs::dir_ls('R'), .f = source)
library(survminer)

dft_ca_ind <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
)
dft_cpt <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_cpt.rds')
)

out_dir <- here('output', 'fig', 'manu')




  



# age_dx is a floored version of ca_cadx_int in this dataset.

lev_age <- c(
  "Age 18-39 at dx.",
  "Age 40-56 at dx."
)

dft_surv_age_bca <- dft_ca_ind %>%
  mutate(
    age_custom = case_when(
      age_dx < 40 & age_dx >= 18 ~ lev_age[1],
      age_dx <= 56 ~ lev_age[2]
    ),
    age_custom = factor(age_custom, levels = lev_age)
  )



dft_surv_age_bca %<>% 
  mutate(
    # in the breast dataset this is fairly simple since stage IV and dmet are
    #   synonymous.
    dx_to_dmets_yrs = if_else(
      stage_dx_iv %in% "Stage IV", 
      0, 
      dx_to_dmets_yrs
    ),
    tt_os_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ tt_os_dx_yrs - dx_to_dmets_yrs,
    )
    # pfs is already relative to stage IV or dmet date.
  ) 

dft_surv_age_bca %<>%
  filter(!is.na(tt_os_dmet_yrs) & !is.na(bca_subtype_f_simple)) %>%
  select(
    record_id, ca_seq, 
    age_custom,
    matches("^bca"),
    dx_to_dmets_yrs,
    tt_os_dmet_yrs,
    os_dmet_status = os_dx_status # death is death - no matter what it's from.
  )



dft_first_cpt <- dft_cpt %>%
  group_by(record_id, ca_seq) %>%
  arrange(cpt_number) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(record_id, ca_seq, dx_cpt_rep_yrs)

dft_surv_age_bca %<>%
  left_join(
    ., dft_first_cpt, by = c('record_id', 'ca_seq')
  )
  


add_specific_dmet_vars <- function(dat) {
  dat %<>% mutate(
    dat,
    tt_cpt_dmet_yrs = case_when(
      is.na(dx_to_dmets_yrs) ~ NA_real_,
      T ~ dx_cpt_rep_yrs - dx_to_dmets_yrs,
    ),
    # glmnet won't take negative times:
    tt_cpt_dmet_yrs_pos = if_else(tt_cpt_dmet_yrs < 0, 0, tt_cpt_dmet_yrs)
  )
  return(dat)
}

dft_surv_age_bca %<>% add_specific_dmet_vars(.)

# ggplot(dft_surv_age_bca, aes(x = tt_cpt_dmet_yrs, y = tt_os_dmet_yrs)) + geom_point()

dft_surv_age_bca %<>%
  filter(tt_cpt_dmet_yrs_pos <= tt_os_dmet_yrs)







# There's a bug in survminer's facet function.  This is a workaround 
#   where we make each plot individually.
surv_age_plot_helper <- function(
    dat,
    bca_subgroup,
    pal = c('#bb5566', '#004488', '#ddaa33')
) {
  dat %<>%
    filter(bca_subtype_f_simple %in% bca_subgroup)
  
  if (nrow(dat) < 1) {
    cli_abort("Filtered down to zero rows with bca_subgroup = {bca_subgroup}")
  }
  n_for_title <- nrow(dat)
  
  surv_obj <- with(
    dat,
    Surv(
      time = tt_cpt_dmet_yrs_pos,
      time2 = tt_os_dmet_yrs,
      event = os_dmet_status
    )
  )
  
  fit <- survfit2(
    surv_obj ~ age_custom,
    data = dat
  )
  
  gg <- ggsurvfit(
    x = fit
  ) + 
    add_risktable(
      risktable_stats = c(
        "n.risk",
        "cum.censor",
        "cum.event"
      ),
      hjust = 0,
      risktable_height = 0.3,
      size = 3.5
    ) +
    add_quantile(
      y_value = 0.5, linetype = 'solid', alpha = 0.3,
      color = "#ddaa33", linewidth = 0.5
    ) + 
    scale_y_continuous(
      expand = c(0,0),
      label = scales::label_percent(),
      name = "Survival"
    ) +
    scale_x_continuous(
      name = "OS from dmet (yrs)",
      expand = expansion(add = 0, mult = c(0, 0.15)), # needed to prevent clipping
      breaks = 0:100
    ) +
    scale_color_manual(
      values = pal
    ) +
    coord_cartesian(
      xlim = c(0, 5.01),
      ylim = c(0,1.01),
      expand = T
    ) +
    labs(
      title = paste0(bca_subgroup, " (n=", n_for_title, ")")
    ) +
    theme(
      axis.title.y = element_blank(),
      plot.title.position = "plot",
      title = element_markdown(),
      # prevents the axis tick label clipping:
      plot.margin=unit(c(.2,.2,.2,.2),"cm")
    )
  
  return(gg)
  
}

gg_hr_pos <- surv_age_plot_helper(dft_surv_age_bca,
                                  bca_subgroup = "HR+, HER2-")
gg_her2_pos <- surv_age_plot_helper(dft_surv_age_bca, 
                                    bca_subgroup = "HER2+")
gg_trip_neg <- surv_age_plot_helper(dft_surv_age_bca, 
                                    bca_subgroup = "Triple Negative")

ggsave(
  plot = gg_hr_pos, height = 4, width = 4,
  filename = here(out_dir, 'fig_age_surv_hr_pos.pdf')
)

ggsave(
  plot = gg_her2_pos, height = 4, width = 4,
  filename = here(out_dir, 'fig_age_surv_her2_pos.pdf')
)

ggsave(
  plot = gg_trip_neg, height = 4, width = 4,
  filename = here(out_dir, 'fig_age_surv_trip_neg.pdf')
)

