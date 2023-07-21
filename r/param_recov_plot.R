plot_param_recov = function(df_chains, df_param){
  p = df_chains %>%
    pivot_longer(cols = b0:b1, names_to = "param", values_to = "value") %>%
    ggplot() +
    stat_halfeye(mapping = aes(x = value, y = param, color = "Posterior")) +
    geom_point(data = df_param,
               mapping = aes(x = sim, y = param, color = "Simulation"), 
               shape = 124, 
               size = 8) +
    scale_color_manual(breaks = c("Simulation", "Posterior"),
                       values = c('Simulation'='red', 'Posterior'='black')) +
    guides(color = guide_legend(title = '',
                                override.aes=list(shape = ""))) +
    theme_minimal() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 15, hjust = 0)) +
    facet_wrap(~param, ncol = 1, scales = "free_x") +
    labs(x = "",
         y = "")
}