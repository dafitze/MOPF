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
source("r/data_simulation.R")
source("r/plot_ce.R")
source("r/plot_pf.R")
source("r/plot_param_recov.R")
d_sim = simulate_data(b0 = 0.0, 
                      b1 = 2.0, 
                      guess = 0.05, 
                      lapse = 0.05, 
                      vpn = 1)
plot_pf(d_sim, mu = 0.0)

d_sim_summary = d_sim %>%
  group_by(stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0)

# model 
# -----------------------------------------------
model = bf(
  response ~ guess + (1 - guess - lapse) * Phi(eta),
  eta ~ 0 + Intercept + stimulus,
  guess ~ 1,
  lapse ~ 1,
  family = bernoulli(link = "identity"),
  nl = TRUE
)

# prior
# ===============================================
# get_prior(model, d_sim)

priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "stimulus", nlpar = "eta"),
  prior(beta(1, 20), nlpar = "lapse", lb = 0, ub = 1),
  prior(beta(1, 20), nlpar = "guess", lb = 0, ub = 1)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = "cmdstanr")

(p_prior_ce = plot_ce(prior_fit, NA, 2))

# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    init = 0,
                    control = list(adapt_delta = 0.99),
                    cores = parallel::detectCores(),
                    backend = "cmdstanr")

(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 2))


# parameter estimates
# -----------------------------------------------
# get_variables(posterior_fit)
chains = posterior_fit %>%
  spread_draws(
    b_eta_Intercept,
    b_eta_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept
  ) %>%
  mutate(
    b0 = b_eta_Intercept,
    b1 = b_eta_stimulus,
    guess = b_guess_Intercept,
    lapse = b_lapse_Intercept
  ) %>%
  select(b0, b1, guess, lapse) 

pars = chains %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi() %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(factor(param, levels = c("b0", "b1", "guess", "lapse"))) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(
           unique(d_sim$b0), 
           unique(d_sim$b1), 
           unique(d_sim$guess), 
           unique(d_sim$lapse)
           )) %>%
  select(param, sim, fit, .lower, .upper)

tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))

# plot overall
# -----------------------------------------------
(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)

