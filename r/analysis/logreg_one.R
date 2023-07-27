library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(ggpmisc)
library(ggpubr)

# ===============================================
# Simulated Data
# ===============================================
source("r/simulation/data_simulation.R")
source("r/plot/plot_ce.R")
source("r/plot/theme_clean.R")
source("r/plot/plot_pf.R")
source("r/plot/plot_priors.R")
source("r/plot/plot_param_recov.R")
source("r/plot/plot_prior_vs_posterior.R")
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
model = bf(response ~ 0 + Intercept + stimulus,
           family = bernoulli('logit'))

# prior
# ===============================================
# get_prior(model, d_sim)
priors = c(
  prior(normal(0.0, 10), class = "b", coef = "stimulus"),
  prior(normal(0.0, 10), class = "b", coef = "Intercept")
)

(p_priors = plot_priors_logreg(priors))

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = 'cmdstanr')

(p_prior_ce = plot_ce(prior_fit, NA, 1))

prior_chains = prior_fit %>%
  spread_draws(
    b_Intercept,
    b_stimulus,
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1 = b_stimulus,
  ) %>%
  select(b0, b1)

# posterior fit
# ===============================================
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 1))

# parameter estimates
# -----------------------------------------------
# get_variables(posterior_fit)
chains = posterior_fit %>%
  spread_draws(
    b_Intercept,
    b_stimulus,
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1 = b_stimulus,
  ) %>%
  select(b0, b1)

pars = chains %>%
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
  mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))


  
# plot overall
# -----------------------------------------------
(p_priors | p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)

plot_prior_vs_posterior(prior_chains, chains)
