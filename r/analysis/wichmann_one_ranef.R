library(brms)
library(cmdstanr)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      b0_sigma = 0.1,
                      b1 = 5.0,
                      b1_sigma = 0.1,
                      b0_guess = 0.1,
                      b0_guess_sigma = 0.01,
                      b0_lapse = 0.1,
                      b0_lapse_sigma = 0.01,
                      n_vpn = 6,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)


# model 
# -----------------------------------------------
model = bf(
  response ~ guess + (1 - guess - lapse) * Phi(eta),
  eta ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
  guess ~ 0 + Intercept + (0 + Intercept | vpn),
  lapse ~ 0 + Intercept + (0 + Intercept | vpn),
  family = bernoulli(link = "identity"),
  nl = TRUE
)

# prior
# ===============================================
# get_prior(model, d_sim)
priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "stimulus", nlpar = "eta"),
  # prior(beta(2, 50), nlpar = "lapse", lb = 0, ub = 1),
  # prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1),
  prior(student_t(3, 0, 2.5), class = "sd", nlpar = "eta"),
  prior(student_t(3, 0, 2.5), class = "b", nlpar = "guess"),
  prior(student_t(3, 0, 2.5), class = "b", nlpar = "lapse"),
  prior(lkj(3), class = "cor")
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = "cmdstanr")

prior_chains = get_chains(prior_fit, type = "wichmann_one_ranef")

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))

# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    init = 0,
                    control = list(adapt_delta = 0.99),
                    cores = parallel::detectCores(),
                    backend = "cmdstanr")

posterior_chains = get_chains(posterior_fit, type = "wichmann_one_ranef")
pars = get_pars(posterior_chains, d_sim)

# tbl = pars %>%
#   ggtexttable(rows = NULL,
#               theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))
# (p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior) | (p_prior_ce / p_posterior_ce))
