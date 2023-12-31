library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulate Data
# ===============================================
d_sim = cell_mean_simulation(b0_pre = 0.0,
                             b1_pre = 2.0,
                             
                             b0_pre_sigma = 0.2,
                             b1_pre_sigma = 0.2,
                             
                             rho = -0.7,
                             
                             n_vpn = 20,
                             n_trials = 20,
                             time = c("pre"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
                             link_function = 'logit')

plot_pf(d_sim, mu = 0.0)


# model
# -----------------------------------------------
model = bf(
  response ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
  family = bernoulli('logit')
)

# prior
# -----------------------------------------------
# get_prior(model, d_sim)
priors = c(
  prior(normal(0, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 5), class = "b", coef = "stimulus"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(3), class = "cor")
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = 'cmdstanr')

prior_chains = prior_fit %>%
  spread_draws(
    b_Intercept,
    b_stimulus,
    sd_vpn__Intercept,
    sd_vpn__stimulus,
  ) |>
  mutate(
    b0_pre = b_Intercept,
    b1_pre = b_stimulus,
    b0_pre_sigma = sd_vpn__Intercept,
    b1_pre_sigma = sd_vpn__stimulus,
  ) |>
  select(b0_pre, b1_pre, b0_pre_sigma, b1_pre_sigma)

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 1, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))


# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

# check fit
# -----------------------------------------------
plot(posterior_fit)
(p_cor = mcmc_pairs(posterior_fit,
           diag_fun = 'dens',
           pars = c('b_Intercept', 'b_stimulus', 'sd_vpn__Intercept', 'sd_vpn__stimulus')))



posterior_chains = posterior_fit %>%
  spread_draws(
    b_Intercept,
    b_stimulus,
    sd_vpn__Intercept,
    sd_vpn__stimulus
  ) |>
  mutate(
    b0_pre = b_Intercept,
    b1_pre = b_stimulus,
    b0_pre_sigma = sd_vpn__Intercept,
    b1_pre_sigma = sd_vpn__stimulus
  ) |>
  select(b0_pre, b1_pre, b0_pre_sigma, b1_pre_sigma)

pars = get_pars(posterior_chains, d_sim)


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 1, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))


# ===============================================
# Repeated Recovery
# ===============================================
rep_recov = 
  repeat_simulation(reps = 50,  # simulate multiple data sets
                    b0_pre = 0.0,
                    b1_pre = 2.0,
                    rho = -0.7,
                    n_vpn = 20,
                    n_trials = 20,
                    time = c("pre"),
                    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |> 
  mutate(
    fit = map(sim_dat, ~update(posterior_fit, newdata = .x)),                    # add fits 
    chains =  map(fit, ~spread_draws(.x, b_Intercept, b_stimulus, sd_vpn__Intercept, sd_vpn__stimulus)),  # add posterior chains
    chains = map(chains, ~select(.x, b0_pre = b_Intercept, b1_pre = b_stimulus, b0_pre_sigma = sd_vpn__Intercept, b1_pre_sigma = sd_vpn__stimulus))) # match names to sim pars
  unnest(pars, names_sep = '_')

(p_rep_recov = plot_rep_recov(rep_recov) + ggtitle("Repetaed Parameter Recovery"))

# plot overall
# -----------------------------------------------
((p_priors | p_posterior)/ p_rep_recov | (p_prior_ce / p_posterior_ce))
