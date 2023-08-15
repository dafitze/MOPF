library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = cell_mean_simulation(b0_pre = 0.0,
                             b0_post = 0.0,
                             
                             b1_pre = 2.0,
                             b1_post = 4.0,
                             
                             n_vpn = 1,
                             n_trials = 20,
                             time = c("pre","post"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)


# model
# -----------------------------------------------
model = bf(
  response ~ inter + predi, 
  lf(inter ~ 0 + time),
  lf(predi ~ 0 + time:stimulus),
  family = bernoulli('probit'),
  nl = TRUE
)

# prior
# ===============================================
get_prior(model, d_sim)
priors = c(
  # Intercepts
  prior(normal(0, 10), nlpar = "inter"),
  
  # Slopes
  prior(gamma(1, 0.1), nlpar = "predi", lb = 0, ub = 5)
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                chains = 4,
                iter = 2000,
                backend = 'cmdstanr')

# prior_chains = get_chains(prior_fit, type = "logreg_prepost")
# get_variables(prior_fit)
prior_chains = prior_fit |>
  spread_draws(
    b_inter_timepre,
    b_inter_timepost,
    `b_predi_timepre:stimulus`,
    `b_predi_timepost:stimulus`
  ) |>
  mutate(
    b0_pre = b_inter_timepre,
    b0_post = b_inter_timepost,
    b1_pre = `b_predi_timepre:stimulus`,
    b1_post = `b_predi_timepost:stimulus`
  ) |>
  select(b0_pre, b0_post, b1_pre, b1_post)

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))

# posterior fit
# ===============================================
posterior_fit = brm(model,
                    data = d_sim,
                    prior = priors,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

# check fit
# -----------------------------------------------
plot(posterior_fit)
# check fit
# -----------------------------------------------
plot(posterior_fit)
(p_cor = mcmc_pairs(posterior_fit,
                    diag_fun = 'dens',
                    pars = c('b_inter_timepre', 'b_inter_timepost', 'b_predi_timepre:stimulus', 'b_predi_timepost:stimulus')))


# posterior_chains = get_chains(posterior_fit, type = "logreg_prepost")
posterior_chains = posterior_fit |>
  spread_draws(
    b_inter_timepre,
    b_inter_timepost,
    `b_predi_timepre:stimulus`,
    `b_predi_timepost:stimulus`
  ) |>
  mutate(
    b0_pre = b_inter_timepre,
    b0_post = b_inter_timepost,
    b1_pre = `b_predi_timepre:stimulus`,
    b1_post = `b_predi_timepost:stimulus`
  ) |>
  select(b0_pre, b0_post, b1_pre, b1_post)


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))

# ===============================================
# Repeated Recovery
# ===============================================
rep_recov = 
  repeat_simulation(reps = 2,  # simulate multiple data sets
                    b0_pre = 0.0,
                    b0_post = 0.0,
                    
                    b1_pre = 2.0,
                    b1_post = 4.0,
                    
                    n_vpn = 1,
                    n_trials = 20,
                    time = c("pre","post"),
                    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |> 
  mutate(
    fit = map(sim_dat, ~update(posterior_fit, newdata = .x)),                    # add fits 
    chains =  map(fit, ~spread_draws(.x, b_timepre
                                     
                                     HIER WEITER
                                     , b_stimulus, sd_vpn__Intercept, sd_vpn__stimulus)),  # add posterior chains
    chains = map(chains, ~select(.x, b0_pre = b_Intercept, b1_pre = b_stimulus, b0_pre_sigma = sd_vpn__Intercept, b1_pre_sigma = sd_vpn__stimulus))) # match names to sim pars
unnest(pars, names_sep = '_')

(p_rep_recov = plot_rep_recov(rep_recov) + ggtitle("Repetaed Parameter Recovery"))

# plot overall
# -----------------------------------------------
((p_priors | p_posterior)/ p_rep_recov | (p_prior_ce / p_posterior_ce))
