plot_lasso_model_comp <- function(
    dat,
    pal = c("#009988", "#ee3377")
    #pal = c("#004488", "#bb5566")
) {
  gg_stab <- ggplot(
    dat,
    aes(x = stability, y = term, color = model)
  ) + 
    geom_point(alpha = 0.8, size = 1.5) + 
    theme_bw() + 
    scale_color_manual(values = pal) + 
    scale_y_discrete(limits = rev) + 
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom"
    ) 
  
  gg_est <- ggplot(
    dat,
    aes(x = log_hr, y = term, color = model)
  ) + 
    geom_point(alpha = 0.8, size = 1.5) + 
    theme_bw() + 
    scale_color_manual(values = pal) + 
    scale_y_discrete(limits = rev) + 
    theme(
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  comb <- cowplot::plot_grid(
    gg_stab, gg_est, nrow = 1, align = 'hv', axis = 'b'
  )
  return(comb)
  
}