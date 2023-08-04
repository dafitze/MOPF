library(brms)
library(cmdstanr)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      b0_sigma = 0.8,
                      b1 = 2.0,
                      b1_sigma = 0.4,
                      rho = -0.7,
                      n_vpn = 6,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

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

prior_chains = get_chains(prior_fit, type = "logreg_one_ranef")

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

posterior_chains = get_chains(posterior_fit, type = "logreg_one_ranef")
pars = get_pars(posterior_chains, d_sim)

# tbl = pars %>%
#   ggtexttable(rows = NULL,
#               theme = ttheme('blank'))




# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 1, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))
# (p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior) | (p_prior_ce / p_posterior_ce))
