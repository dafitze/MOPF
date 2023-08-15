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
                             stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
                             link_function = 'logit')

plot_pf(d_sim, mu = 0.0)

# model
# -----------------------------------------------
model = bf(
  response ~ guess + (1 - guess - lapse) * Phi(inter + preds),
  lf(inter ~ 0 + time),
  lf(preds ~ 0 + time:stimulus),
  # eta ~ 0 + time + time:stimulus,
  guess ~ 0 + Intercept, #time,
  lapse ~ 0 + Intercept, #time,
  nl = TRUE,
  family = bernoulli(link = 'identity')
)

# prior
# ===============================================
priors = c(
  # Priors Psychometric Function
  prior(normal(0, 1), nlpar = "inter"), # Intercepts
  prior(student_t(3 ,0, 10), nlpar = "preds", lb = 0, ub = Inf), # slopes
  # Priors Lapse Rate
  prior(beta(2, 50), class = "b", nlpar = "lapse", lb = 0, ub = 1),
  # Priors Guess Rate
  prior(beta(2, 50), class = "b", nlpar = "guess", lb = 0, ub = 1)
)

prior_fit = brm(model,
                data = d_sim,
                prior = priors,
                sample_prior = "only",
                backend = 'cmdstanr')

# prior_chains = get_chains(prior_fit, type = "wichmann_prepost")
prior_chains = prior_fit %>%
  spread_draws(
    b_inter_timepre,
    b_inter_timepost,
    `b_preds_timepre:stimulus`,
    `b_preds_timepost:stimulus`,
    b_guess_Intercept,
    # b_guess_timepost,
    b_lapse_Intercept,
    # b_lapse_timepost
  ) %>%
  mutate(
    b0_pre = b_inter_timepre,
    b0_post = b_inter_timepost,
    b1_pre = `b_preds_timepre:stimulus`,
    b1_post = `b_preds_timepost:stimulus`,
    guess = b_guess_Intercept,
    # guess_post = b_guess_timepost,
    lapse = b_lapse_Intercept
    # lapse_post = b_lapse_timepost,
  ) %>%
  select(b0_pre, b0_post, b1_pre, b1_post, guess, lapse)# lapse_pre, lapse_post)

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
                    init = 0,
                    control = list(adapt_delta = 0.99), #, max_treedepth = 15),
                    backend = 'cmdstanr')


posterior_chains = posterior_fit %>%
  spread_draws(
    b_inter_timepre,
    b_inter_timepost,
    `b_preds_timepre:stimulus`,
    `b_preds_timepost:stimulus`,
    b_guess_Intercept,
    # b_guess_timepost,
    b_lapse_Intercept,
    # b_lapse_timepost
  ) %>%
  mutate(
    b0_pre = b_inter_timepre,
    b0_post = b_inter_timepost,
    b1_pre = `b_preds_timepre:stimulus`,
    b1_post = `b_preds_timepost:stimulus`,
    guess_pre = b_guess_Intercept,
    # guess_post = b_guess_timepost,
    lapse_pre = b_lapse_Intercept
    # lapse_post = b_lapse_timepost,
  ) %>%
  select(b0_pre, b0_post, b1_pre, b1_post, guess_pre, lapse_pre)# lapse_pre, lapse_post)


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))

# ===============================================
# Repeated Recovery
# ===============================================
rep_recov = 
  repeat_simulation(
    reps = 2,
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
    chains =  map(fit, ~spread_draws(.x, b_inter_timepre, b_inter_timepost, `b_preds_timepre:stimulus`, `b_preds_timepost:stimulus`, b_guess_Intercept, b_lapse_Intercept)),  # add posterior chains
    chains = map(chains, ~select(.x, b0_pre = b_inter_timepre, b0_post = b_inter_timepost, b1_pre = `b_preds_timepre:stimulus`, b1_post = `b_preds_timepost:stimulus`,lapse_pre = b_lapse_Intercept, guess_pre = b_guess_Intercept))) # match names to sim pars


(p_rep_recov = plot_rep_recov(rep_recov) + ggtitle("Repetaed Parameter Recovery"))

# plot overall
# -----------------------------------------------
((p_priors | p_posterior)/ p_rep_recov | (p_prior_ce / p_posterior_ce))
