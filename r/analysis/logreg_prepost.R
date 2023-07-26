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
                              b1 = 2.0, 
                              guess = 0, 
                              lapse = 0, 
                              b2 = 0.0,
                              b3 = 0.8,
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
# model = bf(
#   response ~ 0 + Intercept + time * stimulus,
#   family = bernoulli('logit')
# )
model = bf(
  response ~ 0 + Intercept + time:stimulus,
  family = bernoulli('logit')
)
# prior
# ===============================================
# get_prior(model, d_sim)

priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept"),
  prior(normal(0.0, 10), class = "b", coef = "timepost:stimulus"),
  prior(normal(0.0, 10), class = "b", coef = "timepre:stimulus")
)

# priors = c(
#   prior(normal(0, 10), class = "b", coef = "Intercept"),
#   prior(normal(0.0, 10), class = "b", coef = "stimulus"),
#   prior(normal(0, 1), class = "b", coef = "timepost"),
#   prior(normal(0, 1), class = "b", coef = "timepost:stimulus")
# )

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

(p_prior_ce = plot_ce(prior_fit, NA, 1))

# posterior fit
# ===============================================
posterior_fit = brm(model,
                    data = d_sim,
                    prior = priors,
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

pars = chains %>%
  pivot_longer(cols = 1:4, names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi(.width = 0.93) %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(c("b0", "b1_pre", "b1_post", "b1_diff")) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(unique(d_sim$b0),
                 unique(d_sim$b1_pre),
                 unique(d_sim$b1_post),
                 unique(d_sim$b3)
         )) %>%
  select(param, sim, fit, .lower, .upper, .width)
  

tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

tbl_pars = pars %>%
  # mutate(description = c("pre bias (logit-scale)", 
                         # "pre slope (logit-scale)", 
                         # "additive post bias (logit-scale)", 
                         # "additive post slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))

# plot overall
# -----------------------------------------------
(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl_pars)
