compare_models = function(chains_model_1, chains_model_2, df_param = NA, show.pointin = F){
  pri = chains_model_1 %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(type = "chains 1")
  po = chains_model_2 %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(type = "chains 2")
  
  rbind(pri, po) %>%
    ggplot() +
    stat_halfeye(mapping = aes(x = value, y = 0, fill = type), 
                 .width = c(.51,.93), 
                 alpha = 0.3,
                 # fill = 'orange',
                 # side = 'bottom',
                 geom = "slab") +  
    {if(show.pointin) stat_halfeye(mapping = aes(x = value, y = 0, color = type), .width = c(.51,.93), alpha = 0.3, geom = "pointinterval")} +
    {if(!is.na(df_param)[1])geom_point(data = df_param,
                                       mapping = aes(x = sim, y = 0, color = "Simulation"), 
                                       shape = 124, 
                                       size = 8) } +
    theme_clean() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 15, hjust = 0)) +
    facet_wrap(~param, scales = "free") +
    labs(x = "",
         y = "") +
    # lims(x = c(-20,20)) +
    ggtitle("Prior Chains vs. Posterior Chains")
  
}

plot_prior_vs_posterior = function(prior_chains, posterior_chains){
  pri = prior_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value")
  po = posterior_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value")
  
  ggplot() +
    stat_halfeye(data = pri,
                 mapping = aes(x = value, y = 0),
                 .width = c(.51,.93),
                 alpha = 0.6,
                 fill = 'orange',
                 # side = 'bottom',
                 geom = "slab") +
    stat_halfeye(data = po,
                 mapping = aes(x = value, y = 0),
                 .width = c(.51,.93),
                 alpha = 0.6,
                 fill = 'cyan',
                 geom = "slab") +
    
    theme_clean() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 15, hjust = 0)) +
    facet_wrap(~param, scales = "free") +
    labs(x = "",
         y = "") +
    # lims(x = c(-20,20)) +
    ggtitle("Prior Chains vs. Posterior Chains")
  
  
}