priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "stimulus", nlpar = "eta"),
  prior(beta(1, 20), nlpar = "lapse", lb = 0, ub = 1),
  prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1)
)

plot_priors_wichmann = function(priors){
  p1 = parse_dist(priors) |>
    # filter(.dist == "student_t") |>
    filter(nlpar == "eta") |>
    ggplot(aes(y = 1, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_wrap(~coef, scales = "free", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_tidybayes()
  
  p2 = parse_dist(priors) |>
    # filter(.dist == "beta") |>
    filter(nlpar != "eta") |>
    ggplot(aes(y = 1, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_wrap(~nlpar, scales = "free", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_tidybayes()
  
  p1 | p2
}