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
                             b1_pre = 2.0,
                             
                             lapse_pre = 0.05,
                             guess_pre = 0.05,
                             
                             n_vpn = 1,
                             n_trials = 20,
                             time = c("pre"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)


# model 
# -----------------------------------------------
model = bf(
  response ~ guess + (1 - guess - lapse) * Phi(a + s),
  lf(a ~ 0 + Intercept),
  lf(s ~ 0 + stimulus),
  # eta ~ 0 + Intercept + stimulus,
  guess ~ 1,
  lapse ~ 1,
  family = bernoulli(link = "identity"),
  nl = TRUE
)

# prior
# ===============================================
get_prior(model, d_sim)
priors = c(
  prior(normal(0, 1), class = "b", nlpar = "a"),
  prior(student_t(3 ,0, 2), class = "b", nlpar = "s", lb = 0, ub = Inf),
  prior(beta(2, 50), nlpar = "lapse", lb = 0, ub = 1),
  prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = "cmdstanr")

# prior_chains = get_chains(prior_fit, type = "wichmann_one")
get_variables(prior_fit)
prior_chains = prior_fit %>%
  spread_draws(
    b_a_Intercept,
    b_s_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept
  ) %>%
  mutate(
    b0_pre = b_a_Intercept,
    b1_pre = b_s_stimulus,
    guess_pre = b_guess_Intercept,
    lapse_pre = b_lapse_Intercept
  ) %>%
  select(b0_pre, b1_pre, guess_pre, lapse_pre)

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

# posterior_chains = get_chains(posterior_fit, type = "wichmann_one")
posterior_chains = posterior_fit %>%
  spread_draws(
    b_a_Intercept,
    b_s_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept
  ) %>%
  mutate(
    b0_pre = b_a_Intercept,
    b1_pre = b_s_stimulus,
    guess_pre = b_guess_Intercept,
    lapse_pre = b_lapse_Intercept
  ) %>%
  select(b0_pre, b1_pre, guess_pre, lapse_pre)


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

