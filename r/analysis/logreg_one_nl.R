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
                      b1 = 2.0,
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

# model
# -----------------------------------------------
model = bf(
  response ~ a + s,
  lf(a ~ 0 + Intercept),
  lf(s ~ 0 + stimulus),
  family = bernoulli('logit'),
  nl = TRUE
)

# prior
# ===============================================
priors = c(
  prior(normal(0, 10), nlpar = "a"),
  prior(normal(0, 10), nlpar = "s", lb = 0, ub = Inf)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                chains = 4,
                iter = 20000,
                backend = 'cmdstanr')

# get_variables(prior_fit)
prior_chains = prior_fit %>%
  spread_draws(
    b_a_Intercept,
    b_s_stimulus,
  ) %>%
  mutate(
    b0 = b_a_Intercept,
    b1 = b_s_stimulus,
  ) %>%
  select(b0, b1)

(p_prior_ce = plot_ce(prior_fit, NA, 2))
(p_priors = plot_chains(prior_chains, NA, "orange", "Prior Distributions", F))


# posterior fit
# ===============================================
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

# get_variables(posterior_fit)
posterior_chains = posterior_fit %>%
  spread_draws(
    b_a_Intercept,
    b_s_stimulus,
  ) %>%
  mutate(
    b0 = b_a_Intercept,
    b1 = b_s_stimulus,
  ) %>%
  select(b0, b1)

pars = posterior_chains %>%
  pivot_longer(cols = b0:b1, names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi(.width = 0.93) %>%
  select(param, value, .lower, .upper, .width) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(unique(d_sim$b0),
                 unique(d_sim$b1)
         )) %>%
  select(param, sim, fit, .lower, .upper, .width)

tbl = pars %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 2))
(p_posterior = plot_chains(posterior_chains, pars, "cyan", "Posterior Distributions", T))
(p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))
  
# plot overall
# -----------------------------------------------
((p_priors / p_posterior / p_combo) | (p_prior_ce / p_posterior_ce / tbl))