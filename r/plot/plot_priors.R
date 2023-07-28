# priors = c(
#   prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
#   prior(normal(0, 10), class = "b", coef = "stimulus", nlpar = "eta"),
#   prior(beta(1, 20), nlpar = "lapse", lb = 0, ub = 1),
#   prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1)
# )

# priors = c(
  # prior(normal(0.0, 10), class = "b", coef = "stimulus"),
  # prior(normal(0.0, 10), class = "b", coef = "Intercept")
# )
# priors = c(
#   prior(normal(0.0, 10), class = "b", coef = "Intercept"),
#   prior(normal(0.0, 10), class = "b", coef = "stimulus"),
#   prior(student_t(3, 0, 2.5), class = "sd"),
#   prior(lkj(3), class = "cor")
# )

priors_nl = c(
  prior(normal(0.0, 10), nlpar = "a"),
  prior(normal(0.0, 10), nlpar = "s", lb = 0, ub = Inf)
)

# plot_prior_nl = function(priors){
  parse_dist(priors_nl) |>
    # filter(.dist == "student_t") |>
    # filter(nlpar == "eta") |>
    ggplot(aes(y = class, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_grid(~nlpar, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme_clean() +
    ggtitle("Prior Distributions")
# }

plot_priors_logreg = function(priors){
  parse_dist(priors) |>
    # filter(.dist == "student_t") |>
    # filter(nlpar == "eta") |>
    ggplot(aes(y = class, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_grid(~coef, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme_clean() +
    ggtitle("Prior Distributions")
}

plot_priors_wichmann = function(priors){
  p1 = parse_dist(priors) |>
    # filter(.dist == "student_t") |>
    filter(nlpar == "eta") |>
    ggplot(aes(y = 1, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_wrap(~coef, scales = "free", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_clean() +
    ggtitle("Prior Distributions")
  
  p2 = parse_dist(priors) |>
    # filter(.dist == "beta") |>
    filter(nlpar != "eta") |>
    ggplot(aes(y = 1, dist = .dist, args = .args)) +
    stat_dist_halfeye() +
    facet_wrap(~nlpar, scales = "free", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_clean()
  
  p1 | p2
}