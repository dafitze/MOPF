library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(ggpubr)
library(RMOPF)

# Simulate Data
# ===============================================
d_sim = cell_mean_simulation(b0_pre = 0.0,
                             b1_pre = 2.0,
                             
                             lapse_pre = 0.05,
                             guess_pre = 0.05,
                             
                             n_vpn = 1,
                             n_trials = 20,
                             time = c("pre"),
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
                             link_function = 'logit')

plot_pf(d_sim, mu = 0.0)


# model 
# -----------------------------------------------
model = bf(
  response ~ guess + (1 - guess - lapse) * Phi(inter + preds),
  lf(inter ~ 0 + Intercept),
  lf(preds ~ 0 + stimulus),
  # eta ~ 0 + Intercept + stimulus,
  guess ~ 1,
  lapse ~ 1,
  family = bernoulli(link = "identity"),
  nl = TRUE
)

# prior
# ===============================================
priors = c(
  prior(normal(0, 1), class = "b", nlpar = "inter"),
  prior(student_t(3 ,0, 2), class = "b", nlpar = "preds", lb = 0, ub = Inf),
  prior(beta(2, 50), nlpar = "lapse", lb = 0, ub = 1),
  prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = "cmdstanr")

prior_chains = prior_fit %>%
  spread_draws(
    b_inter_Intercept,
    b_preds_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept
  ) %>%
  mutate(
    b0_pre = b_inter_Intercept,
    b1_pre = b_preds_stimulus,
    guess_pre = b_guess_Intercept,
    lapse_pre = b_lapse_Intercept
  ) %>%
  select(b0_pre, b1_pre, guess_pre, lapse_pre)

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))

# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    init = 0,
                    control = list(adapt_delta = 0.99),
                    cores = parallel::detectCores(),
                    backend = "cmdstanr")

# check fit
# -----------------------------------------------
plot(posterior_fit)
(p_cor = mcmc_scatter(posterior_fit,
                      pars = c('b_inter_Intercept', 'b_preds_stimulus'),
                      size = 1) +
    ggtitle("Correlation Posterior Samples"))

posterior_chains = posterior_fit %>%
  spread_draws(
    b_inter_Intercept,
    b_preds_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept
  ) %>%
  mutate(
    b0_pre = b_inter_Intercept,
    b1_pre = b_preds_stimulus,
    guess_pre = b_guess_Intercept,
    lapse_pre = b_lapse_Intercept
  ) %>%
  select(b0_pre, b1_pre, guess_pre, lapse_pre)


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))

# ===============================================
# Repeated Recovery
# ===============================================
rep_recov = 
  repeat_simulation(
    reps = 50,
    b0_pre = 0.0,
    b1_pre = 2.0,
    lapse_pre = 0.05,
    guess_pre = 0.05,
    n_vpn = 1,
    n_trials = 20,
    time = c("pre"),
    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |> 
  mutate(
    fit = map(sim_dat, ~update(posterior_fit, newdata = .x)),                    # add fits 
    chains =  map(fit, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus)),  # add posterior chains
    chains = map(chains, ~select(.x, b0_pre = b_inter_Intercept, b1_pre = b_preds_stimulus))) # match names to sim pars


(p_rep_recov = plot_rep_recov(rep_recov) + ggtitle("Repetaed Parameter Recovery"))

# plot overall
# -----------------------------------------------
((p_priors | p_posterior)/ p_rep_recov | (p_prior_ce / p_posterior_ce))
