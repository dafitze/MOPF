library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(ggpmisc)
library(ggpubr)
source("r/load_src.R")

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      b1 = 3.5,
                      n_vpn = 1,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

d_sim_summary = d_sim %>%
  group_by(stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0)

pars = tibble(param = c("b0", "b1"),
              sim = c(unique(d_sim$b0), unique(d_sim$b1)))

# model
# -----------------------------------------------
model = bf(
  response ~ 0 + Intercept + stimulus,
  family = bernoulli('logit')
)

model_nl = bf(
  response ~ a + s,
  lf(a ~ 0 + Intercept),
  lf(s ~ 0 + stimulus),
  family = bernoulli('logit'),
  nl = TRUE
)

# prior
# ===============================================
# get_prior(model, d_sim)
priors = c(
  prior(normal(0, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 5), class = "b", coef = "stimulus")
)

priors_nl = c(
  prior(normal(0, 5), nlpar = "a"),
  prior(normal(0, 5), nlpar = "s", lb = 0, ub = Inf)
)

# prior fit
# ===============================================
prior_fit = brm(model, prior = priors, data = d_sim, sample_prior = "only", backend = 'cmdstanr')
prior_chains = prior_fit %>%
  spread_draws(b_Intercept, b_stimulus) %>%
  mutate(b0 = b_Intercept, b1 = b_stimulus) %>%
  select(b0, b1)
(prior_ce = plot_ce(prior_fit, NA, 1))
(p_priors = plot_chains(prior_chains, NA, "orange", "Prior Distributions", F))

prior_fit_nl = brm(model_nl, prior = priors_nl, data = d_sim, sample_prior = "only", backend = 'cmdstanr')
prior_chains_nl = prior_fit_nl %>%
  spread_draws(b_a_Intercept, b_s_stimulus) %>%
  mutate(b0 = b_a_Intercept, b1 = b_s_stimulus) %>%
  select(b0, b1)
(prior_ce_nl = plot_ce(prior_fit_nl, NA, 2))
(p_priors_nl = plot_chains(prior_chains_nl, NA, "orange", "Prior Distributions", F))


# posterior fit
# ===============================================
posterior_fit = brm(model, prior = priors, data = d_sim, cores = parallel::detectCores(), chains = 4, iter = 2000, backend = 'cmdstanr')
posterior_chains = posterior_fit %>%
  spread_draws(b_Intercept, b_stimulus) %>%
  mutate(b0 = b_Intercept, b1 = b_stimulus) %>%
  select(b0, b1)
(posterior_ce = plot_ce(posterior_fit, d_sim_summary, 1))
(p_posterior = plot_chains(posterior_chains, pars, "cyan", "Posterior Distributions", T))

posterior_fit_nl = brm(model_nl, prior = priors_nl, data = d_sim, cores = parallel::detectCores(), chains = 4, iter = 2000, backend = 'cmdstanr')
posterior_chains_nl = posterior_fit_nl %>%
  spread_draws(b_a_Intercept, b_s_stimulus) %>%
  mutate(b0 = b_a_Intercept, b1 = b_s_stimulus) %>%
  select(b0, b1)
(posterior_ce_nl = plot_ce(posterior_fit_nl, d_sim_summary, 2))
(p_posterior_nl = plot_chains(posterior_chains_nl, pars, "cyan", "Posterior Distributions", T))


(prior_combo = compare_models(prior_chains, prior_chains_nl, NA, F))
(posterior_combo = compare_models(posterior_chains, posterior_chains_nl, pars, T))

# plots
# ===============================================
(prior_combo | (prior_ce/prior_ce_nl)) / (posterior_combo | (posterior_ce/posterior_ce_nl))
 
 
 
 

