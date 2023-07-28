plot_chains = function(df_chains, df_param, color, title, show.pointin = F){
  df_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    ggplot() +
    stat_halfeye(mapping = aes(x = value, y = 0), .width = c(.51,.93), fill = color, geom = "slab") +
    {if(show.pointin) stat_halfeye(mapping = aes(x = value, y = 0), .width = c(.51,.93), fill = color, geom = "pointinterval")} +
    {if(!is.na(df_param)[1])geom_point(data = df_param,
                                       mapping = aes(x = sim, y = 0, color = "Simulation"), 
                                       shape = 124, 
                                       size = 8) } +
    facet_wrap(~param, scales = "free") +
    theme_clean() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.title = element_blank(),
          # legend.position = "bottom",
          strip.text.x = element_text(size = 15, hjust = 0)) +
    labs(x = "",
         y = "") +
    ggtitle(title)
}
