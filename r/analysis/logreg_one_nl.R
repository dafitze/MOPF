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
                             b1_pre = 2.5,
                             
                             n_vpn = 1,
                             n_trials = 20,
                             time = c("pre"),
                             stimulus = 
                               # seq(from = -2.14, to = 2.14, length.out = 50)
                               c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)
plot_pf(d_sim, mu = 0.0)

# model
# -----------------------------------------------
model = bf(
  response ~ inter + preds,
  lf(inter ~ 0 + Intercept),
  lf(preds ~ 0 + stimulus),
  family = bernoulli('probit'),
  nl = TRUE
)

# prior
# ===============================================
priors = c(
  prior(normal(0, 1), nlpar = "inter"),
  prior(student_t(3, 0, 2), nlpar = "preds", lb = 0, ub = Inf)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                chains = 4,
                iter = 2000,
                backend = 'cmdstanr')


# prior_chains = get_chains(prior_fit, type = "logreg_one_nl")
prior_chains = prior_fit |>
  spread_draws(
    b_inter_Intercept,
    b_preds_stimulus
  ) |>
  mutate(
    b0_pre = b_inter_Intercept,
    b1_pre = b_preds_stimulus
  ) |>
  select(b0_pre, b1_pre)

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))


# posterior fit
# ===============================================
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

# check fit
# ===============================================
plot(posterior_fit)
# mcmc_pairs(posterior_fit,
#            diag_fun = 'dens',
#            np = nuts_params(posterior_fit))
(p_cor = mcmc_scatter(posterior_fit,
             pars = c('b_inter_Intercept', 'b_preds_stimulus'),
             size = 1) +
    ggtitle("Correlation Posterior Samples"))

brms::rhat(posterior_fit) |>
  mcmc_rhat() + 
  yaxis_text(hjust = 1)



# posterior_chains = get_chains(posterior_fit, type = "logreg_one_nl")
posterior_chains = posterior_fit |>
  spread_draws(
    b_inter_Intercept,
    b_preds_stimulus
  ) |>
  mutate(
    b0_pre = b_inter_Intercept,
    b1_pre = b_preds_stimulus
  ) |>
  select(b0_pre, b1_pre)

pars = get_pars(posterior_chains, d_sim)

# tbl = pars %>%
#   ggtexttable(rows = NULL,
#               theme = ttheme('blank'))


# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))
(p_combo = plot_chains(list(prior = prior_chains, posterior = posterior_chains),
                       plot_data = d_sim,
                       title = "Prior vs. Posterior"))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior / p_cor) | (p_prior_ce / p_posterior_ce))

# ===============================================
# Repeated Recovery
# ===============================================
rep_recov = 
  # simulate multiple data sets
  repeat_simulation(reps = 20,
                    
                    b0_pre = 0.0,
                    b1_pre = 2.0,
                    
                    n_vpn = 1,
                    n_trials = 20,
                    time = c("pre"),
                    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |> 
  mutate(
    # add fits  
    fit = map(sim_dat, ~update(posterior_fit, newdata = .x)),
    # add posterior chains
    chains =  map(fit, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus)),
    # rename chains to match simulation pars names
    chains = map(chains, ~select(.x, b0_pre = b_inter_Intercept, b1_pre = b_preds_stimulus)),
    # add parameters
    pars = map2(chains, sim_dat, ~get_pars(.x, .y))
  ) |>
  unnest(pars, names_sep = '_')


# Plot Repeated Recovery
# -----------------------------------------------
(p_rep_recov = plot_rep_recov(rep_recov))
