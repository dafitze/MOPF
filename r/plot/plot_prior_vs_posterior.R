plot_prior_vs_posterior = function(prior_chains, posterior_chains){
  pri = prior_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value")
  posterior_chains = chains
  po = posterior_chains %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value")
  
  ggplot() +
    stat_halfeye(data = pri,
                 mapping = aes(x = value, y = -0.05), 
                 .width = c(.51,.93), 
                 alpha = 0.6, 
                 fill = 'orange',
                 side = 'bottom') +
    stat_halfeye(data = po,
                 mapping = aes(x = value, y = 0.05), .width = c(.51,.93), alpha = 0.6, fill = 'cyan') +
    # stat_halfeye( mapping = aes(x = value, y = 1.5, color = "Prior"), .width = c(.51,.93)) +
    # geom_point(data = pars,
    # mapping = aes(x = sim, y = 1.5, color = "Simulation"), 
    # shape = 124, 
    # size = 8) +
    # scale_color_manual(breaks = c("Simulation", "Posterior"),
    # values = c('Simulation'='red', 'Posterior'='black', 'Prior'='orange')) +
    # guides(color = guide_legend(title = '',
    # override.aes=list(shape = ""))) +
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