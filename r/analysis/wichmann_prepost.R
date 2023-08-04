library(brms)
library(cmdstanr)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      b1 = 2.0,
                      b2 = 0.0,
                      b3 = 2.0,
                      b0_guess = 0.1,
                      b1_guess = 0.01,
                      b0_lapse = 0.1,
                      b1_lapse = 0.01,
                      n_vpn = 1,
                      n_trials = 20,
                      time = c("pre","post"),
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

# model
# -----------------------------------------------
model = bf(
  response ~ Phi(guess) + Phi(1 - guess - lapse) * Phi(eta),
  eta ~ 0 + Intercept + time:stimulus,
  guess ~ 0 + Intercept + time,
  lapse ~ 0 + Intercept + time,
  nl = TRUE,
  family = bernoulli(link = 'identity')
)

# prior
# ===============================================
# get_prior(model, d_sim)
priors = c(
  # Priors Psychometric Function
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "timepre:stimulus", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "timepost:stimulus", nlpar = "eta"),
  
  # Priors Lapse Rate
  prior(beta(2, 50), nlpar = "lapse", lb = 0, ub = 1),
  # prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "lapse"),
  prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "lapse"),
  
  # Priors Guess Rate
  prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1),
  # prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "guess"),
  prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "guess")
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

prior_chains = get_chains(prior_fit, type = "wichmann_prepost")


(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))

# posterior fit
# ===============================================
posterior_fit = brm(model,
                    data = d_sim,
                    prior = priors,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    init = 0,
                    control = list(adapt_delta = 0.99), #, max_treedepth = 15),
                    backend = 'cmdstanr')

posterior_chains = get_chains(posterior_fit, type = "wichmann_prepost")
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

