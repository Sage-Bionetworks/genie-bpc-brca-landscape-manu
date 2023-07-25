plot_hr_rel_scatter <- function(
    dat,
    legend_title = NULL,
    plot_subtitle = NULL,
    plot_title = NULL,
    pal = c("#ddaa33", "#bb5566", "#004488", "gray80"),
    txt_size = 2.5) {
  
  dat %<>%
    mutate(
      type = case_when(
        str_detect(feature, "_mut$") ~ "Mutation",
        str_detect(feature, "_cna$") ~ "CNA",
        str_detect(feature, "_fus$") ~ "Fusion",
        T ~ "Other"
      ),
      type = forcats::fct_inorder(type),
      feature = as.character(feature),
      feature = str_replace(feature, 
                            pattern = "_mut|_cna|_fus",
                            replacement = "")
    )
  
  gg <- ggplot(data = dat,
               aes(x = hr, y = reliability, color = type)) +
    geom_vline(xintercept = 1, linewidth = 2, color = "#009988", alpha = 0.2) + 
    geom_point() +
    geom_text_repel(aes(label = feature), color = "black",
                    size = txt_size) + 
    theme_bw() +
    scale_color_manual(name = legend_title,
                       values = pal) + 
    labs(title = plot_title,
         subtitle = plot_subtitle, 
         x = "Cumulative hazard ratio",
         y = "Reliability") + 
    scale_y_continuous(limits = c(NA, 1)) + 
    theme(
      axis.title.y = element_markdown(vjust = 0.5, angle = 0),
      plot.title.position = "plot"
    )
  
  
  return(gg)
  
}