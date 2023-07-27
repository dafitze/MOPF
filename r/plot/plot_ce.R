plot_ce = function(fit, df_summary = NA, index = 1){
  if (is_tibble(df_summary)){
    plot(conditional_effects(fit), plot = F)[[index]] +
      geom_point(data = df_summary, aes(stimulus, mean_response), size = 3, alpha = 1.0) +
      lims(x = c(-2,2), y = c(0,1)) +
      theme_clean() +
      ggtitle("posterior predictive")
  }
  else {
    plot(conditional_effects(fit), plot = F)[[index]] +
      lims(x = c(-2,2), y = c(0,1)) +
      theme_clean() +
      ggtitle("prior predictive")
  }
  
}