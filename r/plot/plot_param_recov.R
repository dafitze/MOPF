plot_param_recov = function(df_chains, df_param){
  p = df_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    # arrange(c("b0", "b1_pre", "b1_post", "b1_diff")) %>%
    ggplot() +
    stat_halfeye(mapping = aes(x = value, y = 1.5, color = "Posterior"), .width = c(.51,.93)) +
    geom_point(data = df_param,
               mapping = aes(x = sim, y = 1.5, color = "Simulation"), 
               shape = 124, 
               size = 8) +
    scale_color_manual(breaks = c("Simulation", "Posterior"),
                       values = c('Simulation'='red', 'Posterior'='black')) +
    guides(color = guide_legend(title = '',
                                override.aes=list(shape = ""))) +
    theme_clean() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 15, hjust = 0)) +
    facet_wrap(~param, scales = "free") +
    labs(x = "",
         y = "") +
    ggtitle("Parameter Recovery")
}