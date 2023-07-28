library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(ggpmisc)
library(ggpubr)
source("r/load_src.R")

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      # b0_sigma = 0.1,
                      b1 = 5.0,
                      # b1_sigma = 0.1,
                      b0_guess = 0.1,
                      # b0_guess_sigma = 0.01,
                      b0_lapse = 0.1,
                      # b0_lapse_sigma = 0.01,
                      n_vpn = 6,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

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
  eta ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
  guess ~ 0 + Intercept + (0 + Intercept | vpn),
  lapse ~ 0 + Intercept + (0 + Intercept | vpn),
  family = bernoulli(link = "identity"),
  nl = TRUE
)

# prior
# ===============================================
# get_prior(model, d_sim)

priors = c(
  prior(normal(0, 10), class = "b", coef = "Intercept", nlpar = "eta"),
  prior(normal(0, 10), class = "b", coef = "stimulus", nlpar = "eta"),
  prior(beta(2, 50), nlpar = "lapse", lb = 0, ub = 1),
  prior(beta(2, 50), nlpar = "guess", lb = 0, ub = 1),
  prior(student_t(3, 0, 2.5), class = "sd", nlpar = "eta"),
  # prior(student_t(3, 0, 2.5), class = "sd", nlpar = "guess"),
  # prior(student_t(3, 0, 2.5), class = "sd", nlpar = "lapse"),
  prior(lkj(3), class = "cor")
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = "cmdstanr")

prior_chains = prior_fit %>%
  spread_draws(
    b_eta_Intercept,
    b_eta_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept,
    sd_vpn__eta_Intercept,
    sd_vpn__eta_stimulus,
    sd_vpn__guess_Intercept,
    sd_vpn__lapse_Intercept,
    cor_vpn__eta_Intercept__eta_stimulus
  ) %>%
  mutate(
    b0 = b_eta_Intercept,
    b0_sigma = sd_vpn__eta_Intercept,
    b1 = b_eta_stimulus,
    b1_sigma = sd_vpn__eta_stimulus,
    b0_guess = b_guess_Intercept,
    b0_guess_sigma = sd_vpn__guess_Intercept,
    b0_lapse = b_lapse_Intercept,
    b0_lapse_sigma = sd_vpn__lapse_Intercept,
    rho = cor_vpn__eta_Intercept__eta_stimulus
  ) %>%
  select(b0, b0_sigma, b1, b1_sigma, b0_guess, b0_guess_sigma, b0_lapse, b0_lapse_sigma, rho) 


(p_prior_ce = plot_ce(prior_fit, NA, 2))
(p_priors = plot_chains(prior_chains, NA, "orange", "Prior Distributions", F))

# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    init = 0,
                    control = list(adapt_delta = 0.99),
                    cores = parallel::detectCores(),
                    backend = "cmdstanr")

posterior_chains = posterior_fit %>%
  spread_draws(
    b_eta_Intercept,
    b_eta_stimulus,
    b_guess_Intercept,
    b_lapse_Intercept,
    sd_vpn__eta_Intercept,
    sd_vpn__eta_stimulus,
    sd_vpn__guess_Intercept,
    sd_vpn__lapse_Intercept,
    cor_vpn__eta_Intercept__eta_stimulus
  ) %>%
  mutate(
    b0 = b_eta_Intercept,
    b0_sigma = sd_vpn__eta_Intercept,
    b1 = b_eta_stimulus,
    b1_sigma = sd_vpn__eta_stimulus,
    b0_guess = b_guess_Intercept,
    b0_guess_sigma = sd_vpn__guess_Intercept,
    b0_lapse = b_lapse_Intercept,
    b0_lapse_sigma = sd_vpn__lapse_Intercept,
    rho = cor_vpn__eta_Intercept__eta_stimulus
  ) %>%
  select(b0, b0_sigma, b1, b1_sigma, b0_guess, b0_guess_sigma, b0_lapse, b0_lapse_sigma, rho) 

pars = posterior_chains %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  mean_qi(.width = .93) %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(factor(param, levels = c("b0", "b0_sigma", "b1", "b1_sigma", "b0_guess", "b0_guess_sigma", "b0_lapse", "b0_lapse_sigma"))) %>%
  mutate(fit = round(value, digits = 3),
         .lower = round(.lower, digits = 3),
         .upper = round(.upper, digits = 3),
         sim = c(
           unique(d_sim$b0), 
           unique(d_sim$b0_sigma),
           unique(d_sim$b1), 
           unique(d_sim$b1_sigma),
           unique(d_sim$b0_guess),
           unique(d_sim$b0_guess_sigma),
           unique(d_sim$b0_lapse),
           unique(d_sim$b0_lapse_sigma),
           unique(d_sim$rho)
           )) %>%
  select(param, sim, fit, .lower, .upper)

tbl = pars %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 2))
(p_posterior = plot_chains(posterior_chains, pars, "cyan", "Posterior Distributions", T))
(p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior / p_combo) | (p_prior_ce / p_posterior_ce / tbl))
