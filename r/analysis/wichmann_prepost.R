library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = cell_mean_simulation(b0_pre = 0.0,
                             b0_post = 0.0,
                             
                             b1_pre = 2.0,
                             b1_post = 4.0,
                             
                             n_vpn = 1,
                             n_trials = 20,
                             time = c("pre","post"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

# model
# -----------------------------------------------
model = bf(
  response ~ Phi(guess) + Phi(1 - guess - lapse) * Phi(eta),
  eta ~ 0 + time + time:stimulus,
  guess ~ 0 + time,
  lapse ~ 0 + time,
  nl = TRUE,
  family = bernoulli(link = 'identity')
)

# prior
# ===============================================
get_prior(model, d_sim)
priors = c(
  # Priors Psychometric Function
  prior(normal(0, 10), class = "b", coef = "timepre", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "timepost", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "timepre:stimulus", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "timepost:stimulus", nlpar = "eta"),
  
  # Priors Lapse Rate
  prior(beta(2, 50), class = "b", nlpar = "lapse", lb = 0, ub = 1),
  # prior(beta(2, 50), class = "b", coef = "timepost", nlpar = "lapse", lb = 0, ub = 1),
  # prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "lapse"),
  # prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "lapse"),
  
  # Priors Guess Rate
  prior(beta(2, 50), class = "b", nlpar = "guess", lb = 0, ub = 1)
  # prior(beta(2, 50), class = "b", coef = "timepost", nlpar = "guess", lb = 0, ub = 1)
  # prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "guess"),
  # prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "guess")
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

# prior_chains = get_chains(prior_fit, type = "wichmann_prepost")
prior_chains = prior_fit %>%
  spread_draws(
    b_eta_timepre,
    b_eta_timepost,
    `b_eta_timepre:stimulus`,
    `b_eta_timepost:stimulus`,
    b_guess_timepre,
    b_guess_timepost,
    b_lapse_timepre,
    b_lapse_timepost
  ) %>%
  mutate(
    b0_pre = b_eta_timepre,
    b0_post = b_eta_timepost,
    b1_pre = `b_eta_timepre:stimulus`,
    b1_post = `b_eta_timepost:stimulus`,
    guess_pre = b_guess_timepre,
    guess_post = b_guess_timepost,
    lapse_pre = b_lapse_timepre,
    lapse_post = b_lapse_timepost,
  ) %>%
  select(b0_pre, b0_post, b1_pre, b1_post, guess_pre, guess_post, lapse_pre, lapse_post)

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

# posterior_chains = get_chains(posterior_fit, type = "wichmann_prepost")
posterior_chains = posterior_fit %>%
  spread_draws(
    b_eta_timepre,
    b_eta_timepost,
    `b_eta_timepre:stimulus`,
    `b_eta_timepost:stimulus`,
    b_guess_timepre,
    b_guess_timepost,
    b_lapse_timepre,
    b_lapse_timepost
  ) %>%
  mutate(
    b0_pre = b_eta_timepre,
    b0_post = b_eta_timepost,
    b1_pre = `b_eta_timepre:stimulus`,
    b1_post = `b_eta_timepost:stimulus`,
    guess_pre = b_guess_timepre,
    guess_post = b_guess_timepost,
    lapse_pre = b_lapse_timepre,
    lapse_post = b_lapse_timepost,
  ) %>%
  select(b0_pre, b0_post, b1_pre, b1_post, guess_pre, guess_post, lapse_pre, lapse_post)


pars = get_pars(posterior_chains, d_sim)

# tbl = pars %>%
#   ggtexttable(rows = NULL,
#               theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))
(p_combo = plot_chains(list(prior = prior_chains, posterior = posterior_chains),
                       title = "Prior vs. Posterior"))


# plot overall
# -----------------------------------------------
((p_priors / p_posterior) | (p_prior_ce / p_posterior_ce))

