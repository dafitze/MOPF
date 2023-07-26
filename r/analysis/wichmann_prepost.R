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
d_sim = simulate_prepost_data(b0 = 0.0, 
                              b1 = 5.0, 
                              b2 = 0.0,
                              b3 = 0.0,
                              b0_guess = 0.0, 
                              b2_guess = 0.25,
                              b0_lapse = 0.0, 
                              b2_lapse = 0.25,
                              vpn = 1,
                              time = c("pre", "post"))
plot_pf_prepost(d_sim, mu = 0.0)

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
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0.0, 10), class = "b", coef = "timepost:stimulus", nlpar = "eta"),
  prior(normal(0.0, 10), class = "b", coef = "timepre:stimulus", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "lapse"),
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "guess"),
  prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "lapse"),
  prior(normal(0, 1), class = "b", coef = "timepost", nlpar = "guess")
)

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
get_variables(posterior_fit)
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

pars = chains %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi() %>%
  select(param, value, .lower, .upper) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(
           unique(d_sim$b0),
           unique(d_sim$b3),
           unique(d_sim$b1_post),
           unique(d_sim$b1_pre),
           NA, 
           NA,
           NA,
           NA,
           NA,
           NA
         )) %>%
  select(param, sim, fit, .lower, .upper)


tbl = pars %>%
  # mutate(description = c(
  #   "pre intercept (psychometric function)", 
  #   "pre intercept (guess)",
  #   "pre intercept (lapse)",
  #   "pre slope (psychometric function)", 
  #   "additive post intercept (psychometric function)",
  #   "additive post intercept (guess)",
  #   "additive post intercept (lapse)",
  #   "additive post slope (psychometric function)"
  # )) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))

# plot overall
# -----------------------------------------------
(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)






# ===============================================
# PL Data
# ===============================================
# - what changes with the roll stimuli?

# data
# -----------------------------------------------
d = read_csv("../data/exp_pro/d_test.csv") %>%
  mutate(axis = as_factor(axis),
         time = as_factor(measurement),
         vpn = as_factor(vpn),
         trained = as_factor(axis_trained_bool)) %>%
  select(vpn, time, axis, stimulus = stim, response = resp, trained) %>%
  filter(vpn == 3,
         time != "mid", 
         axis == 1)

# summarised vp data
vp_data = d %>%
  group_by(time, stimulus) %>%
  summarise(mean = mean(response)) %>%
  mutate(lower__ = 0, # unused but necessary
         upper__ = 0,
         effect2__ = 0) # unsued but necessary

# posterior fit
# ===============================================
pl_fit = brm(model,
             data = d,
             prior = priors,
             cores = parallel::detectCores(),
             chains = 4,
             iter = 2000,
             init = 0,
             control = list(adapt_delta = 0.99), #, max_treedepth = 15),
             backend = 'cmdstanr')

pl_posterior_ce = plot_ce(pl_fit, d_sim_summary, 4)

# parameter estimates
# -----------------------------------------------
# get_variables(posterior_fit)
pl_pars = get_pars_mixture_prepost(pl_fit)

pl_tbl = pl_pars %>%
  mutate(description = c(
    "pre intercept (psychometric function)", 
    "pre intercept (guess)",
    "pre intercept (lapse)",
    "pre slope (psychometric function)", 
    "additive post intercept (psychometric function)",
    "additive post intercept (guess)",
    "additive post intercept (lapse)",
    "additive post slope (psychometric function)"
  )) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# plot overall
# -----------------------------------------------
(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)
