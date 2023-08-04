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
                      b2 = 0.0,
                      b3 = 2.0,
                      n_vpn = 1,
                      n_trials = 20,
                      time = c("pre","post"),
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

d_sim_summary = d_sim %>%
  group_by(time, stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0,
         effect2__ = time)

# model
# -----------------------------------------------
model = bf(
  response ~ 0 + Intercept + time:stimulus,
  family = bernoulli('logit')
)

# prior
# ===============================================
get_prior(model, d_sim)
priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept"),
  prior(normal(0.0, 10), class = "b", coef = "timepost:stimulus"),
  prior(normal(0.0, 10), class = "b", coef = "timepre:stimulus")
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                chains = 4,
                iter = 2000,
                backend = 'cmdstanr')
summary(prior_fit)

prior_chains = prior_fit %>%
  spread_draws(
    b_Intercept,
    `b_timepre:stimulus`,
    `b_timepost:stimulus`
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1_pre = `b_timepre:stimulus`,
    b1_post = `b_timepost:stimulus`,
    b1_diff = b1_post - b1_pre
  ) %>%
  select(b0, b1_pre, b1_post, b1_diff)

(p_prior_ce = plot_ce(prior_fit, NA, 1))
(p_priors = plot_chains(prior_chains, NA, "orange", "Prior Distributions", F))

# posterior fit
# ===============================================
posterior_fit = brm(model,
                    data = d_sim,
                    prior = priors,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

# get_variables(posterior_fit)
posterior_chains = posterior_fit %>%
  spread_draws(
    b_Intercept,
    `b_timepre:stimulus`,
    `b_timepost:stimulus`
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1_pre = `b_timepre:stimulus`,
    b1_post = `b_timepost:stimulus`,
    b1_diff = b1_post - b1_pre
  ) %>%
  select(b0, b1_pre, b1_post, b1_diff)

pars = posterior_chains %>%
  pivot_longer(cols = 1:4, names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi(.width = 0.93) %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(c("b0", "b1_pre", "b1_post", "b1_diff")) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(unique(d_sim$b0),
                 unique(d_sim$b1),
                 unique(d_sim$b1) + unique(d_sim$b3),
                 unique(d_sim$b3)
         )) %>%
  select(param, sim, fit, .lower, .upper, .width)
  
tbl = pars %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 1))
(p_posterior = plot_chains(posterior_chains, pars, "cyan", "Posterior Distributions", T))
(p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior / p_combo) | (p_prior_ce / p_posterior_ce / tbl))
