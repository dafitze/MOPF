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
                             
                             n_vpn = 1,
                             n_trials = 800,
                             time = c("pre"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

# d_sim = simulate_data(b0 = 0.0,
#                       b1 = 3.5,
#                       n_vpn = 1,
#                       n_trials = 20,
#                       time = "pre",
#                       stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

# model
# -----------------------------------------------
model = bf(
  response ~ 0 + Intercept + stimulus,
  family = bernoulli('logit') # probit üblicher
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
get_prior(model_nl, d_sim)
priors = c(
  prior(normal(0, 1), class = "b", coef = "Intercept"),
  prior(student_t(3, 0, 2), class = "b", coef = "stimulus")
)

priors_nl = c(
  prior(normal(0, 1), nlpar = "a"),
  prior(student_t(3, 0, 2), nlpar = "s", lb = 0, ub = Inf)
)

# prior fit
# ===============================================
prior_fit = brm(model, prior = priors, data = d_sim, sample_prior = "only", backend = 'cmdstanr')
prior_chains = prior_fit %>%
  spread_draws(b_Intercept, b_stimulus) %>%
  mutate(b0_pre = b_Intercept, b1_pre = b_stimulus) %>%
  select(b0_pre, b1_pre)
(prior_ce = plot_ce(prior_fit, plot_data = NA, index = 1, title = "Prior CE linear"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))

prior_fit_nl = brm(model_nl, prior = priors_nl, data = d_sim, sample_prior = "only", backend = 'cmdstanr')
prior_chains_nl = prior_fit_nl %>%
  spread_draws(b_a_Intercept, b_s_stimulus) %>%
  mutate(b0_pre = b_a_Intercept, b1_pre = b_s_stimulus) %>%
  select(b0_pre, b1_pre)
(prior_ce_nl = plot_ce(prior_fit_nl, plot_data = NA, index = 2, title = "Prior CE non-linear"))
(p_priors_nl = plot_chains(prior_chains_nl, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))



# posterior fit
# ===============================================
posterior_fit = brm(model, prior = priors, data = d_sim, cores = parallel::detectCores(), chains = 4, iter = 2000, backend = 'cmdstanr')
posterior_chains = posterior_fit %>%
  spread_draws(b_Intercept, b_stimulus) %>%
  mutate(b0_pre = b_Intercept, b1_pre = b_stimulus) %>%
  select(b0_pre, b1_pre)
(posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 1, title = "Posterior CE linear"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = "cyan", title ="Posterior Distributions", show_pointinterval = T))

posterior_fit_nl = brm(model_nl, prior = priors_nl, data = d_sim, cores = parallel::detectCores(), chains = 4, iter = 2000, backend = 'cmdstanr')
posterior_chains_nl = posterior_fit_nl %>%
  spread_draws(b_a_Intercept, b_s_stimulus) %>%
  mutate(b0_pre = b_a_Intercept, b1_pre = b_s_stimulus) %>%
  select(b0_pre, b1_pre)
(posterior_ce_nl = plot_ce(posterior_fit_nl, plot_data = d_sim, index = 2, title = "Posterior CE non-linear"))
(p_posterior_nl = plot_chains(posterior_chains_nl, plot_data = d_sim, color = "cyan", title = "Posterior Distributions", show_pointinterval = T))


(prior_combo = plot_chains(list(linear = prior_chains, non_linear = prior_chains_nl), 
                           plot_data = NA, 
                           title = "Prior Comparison",
                           show_pointinterval = F))

(posterior_combo = plot_chains(list(linear = posterior_chains, non_linear = posterior_chains_nl), 
                                  plot_data = d_sim, 
                                  title = "Posterior Comparison",
                                  show_pointinterval = T))



# plots
# ===============================================
(prior_combo | (prior_ce/prior_ce_nl)) / (posterior_combo | (posterior_ce/posterior_ce_nl))
 
 
 
 

