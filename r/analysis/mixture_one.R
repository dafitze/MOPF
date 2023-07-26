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
source("02_R/psychometric_function/data_simulation.R")
d_sim = simulate_prepost_data(mu = 0.0, 
                              sigma = 50, 
                              guess = 0, 
                              lapse = 0, 
                              b2 = 0.0,
                              b3 = 0.01,
                              vpn = 1, 
                              time = c("pre", "post"))
plot_pf_prepost(d_sim, mu = 0.5)

d_sim_summary = d_sim %>%
  group_by(stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0)



# model
# -----------------------------------------------
mix = mixture(bernoulli(link = "logit"), 
              bernoulli(link = "logit"),
              order = 'none')
model = bf(
  response ~ 1,                   # initializing dependent variable 
  mu1 ~ 0 + Intercept + stimulus, # psychometric function
  mu2 ~ 0 + Intercept,            # guessing probability
  theta1 ~ 0 + Intercept,         # mixing proportion (lapse rate)
  family = mix
)

# prior
# ===============================================
get_prior(model, d_sim)

priors = c(
  prior(student_t(3, 0, 2.5), class = "b", coef = "Intercept", dpar = "mu1"),
  prior(student_t(3, 0, 10), class = "b", coef = "stimulus", dpar = "mu1"),
  prior(student_t(3, 0, 2.5), class = "b", coef = "Intercept", dpar = "mu2"),
  prior(logistic(10,100), class = "b", coef = "Intercept", dpar = "theta1")
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

p_prior_ce = plot_ce(prior_fit, NA, 1)

# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                  data = d_sim,
                  prior = priors,
                  cores = parallel::detectCores(),
                  chains = 4,
                  iter = 2000,
                  init = 0,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  backend = 'cmdstanr',
                  file = "mixture_one")

p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 1)

# parameter estimates
# -----------------------------------------------
pars = get_pars_mixture_one(posterior_fit)

tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
p_param_recov = plot_param_recov(pars) 

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
         time == "pre", 
         axis == 1)

# summarised vp data
vp_data = d %>%
  group_by(time, stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0, # unused but necessary
         upper__ = 0) # unsued but necessary

# posterior fit
# ===============================================
pl_fit = brm(model,
             data = d,
             prior = priors,
             cores = parallel::detectCores(),
             chains = 4,
             iter = 2000,
             init = 0,
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             backend = 'cmdstanr')

pl_posterior_ce = plot_ce(pl_fit, vp_data, 1)

# parameter estimates
# -----------------------------------------------
pl_pars = get_pars_mixture_one(pl_fit)

pl_tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# plot overall
# -----------------------------------------------
(p_prior_ce | pl_posterior_ce) / pl_tbl
p_posterior_ce | pl_posterior_ce
