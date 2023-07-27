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
source("r/plot/theme_clean.R")
source("r/plot/plot_ce.R")
source("r/plot/plot_pf.R")
source("r/plot/plot_priors.R")
source("r/plot/plot_param_recov.R")
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

d_sim_summary = d_sim %>%
  group_by(time, stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0,
         effect2__ = time)

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

plot_priors_wichmann(priors)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

(p_prior_ce = plot_ce(prior_fit, NA, 2))

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

(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 2))

# parameter estimates
# -----------------------------------------------
# get_variables(posterior_fit)
chains = posterior_fit %>%
  spread_draws(
    b_eta_Intercept,
    `b_eta_timepre:stimulus`,
    `b_eta_timepost:stimulus`,
    b_guess_Intercept,
    b_guess_timepost,
    b_lapse_Intercept,
    b_lapse_timepost
  ) %>%
  mutate(
    b0 = b_eta_Intercept,
    b1_pre = `b_eta_timepre:stimulus`,
    b1_post = `b_eta_timepost:stimulus`,
    b3 = `b_eta_timepost:stimulus`,
    b0_guess = b_guess_Intercept,
    b2_guess = b_guess_timepost,
    b0_lapse = b_lapse_Intercept,
    b2_lapse = b_lapse_timepost,
    b1_diff = b1_post - b1_pre,
    lapse_pre = b0_lapse,
    lapse_post = b0_lapse + b2_lapse,
    lapse_diff = lapse_post - lapse_pre,
    guess_pre = b0_guess,
    guess_post = b0_guess + b2_guess,
    guess_diff = guess_post - guess_pre
  ) %>%
  select(b0, b1_pre, b1_post, b1_diff, guess_pre, guess_post, guess_diff, lapse_pre, lapse_post, lapse_diff)

chains %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi() %>%
  select(param, value, .lower, .upper) %>%
  arrange(factor(param, levels = c("b0", "b1_pre", "b1_post", "b1_diff", "guess_pre", "guess_post", "guess_diff", "lapse_pre", "lapse_post", "lapse_diff")))
  
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(
           unique(d_sim$b0),
           unique(d_sim$b1),
           unique(d_sim$b1) + unique(d_sim$b3),
           unique(d_sim$b3),
           unique(d_sim$b0_guess), 
           unique(d_sim$b0_guess) + unique(d_sim$b1_guess),
           unique(d_sim$b1_guess),
           unique(d_sim$b0_lapse), 
           unique(d_sim$b0_lapse) + unique(d_sim$b1_lapse),
           unique(d_sim$b1_lapse)
         )) %>%
  select(param, sim, fit, .lower, .upper)


tbl = pars %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))

# plot overall
# -----------------------------------------------
(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)

