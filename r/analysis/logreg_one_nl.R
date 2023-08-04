library(brms)
library(cmdstanr)
library(patchwork)
library(ggpubr)
library(RMOPF)

# ===============================================
# Simulated Data
# ===============================================
d_sim = simulate_data(b0 = 0.0,
                      b1 = 2.0,
                      n_vpn = 1,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)


# model
# -----------------------------------------------
model = bf(
  response ~ a + s,
  lf(a ~ 0 + Intercept),
  lf(s ~ 0 + stimulus),
  family = bernoulli('logit'),
  nl = TRUE
)

# prior
# ===============================================
priors = c(
  prior(normal(0, 10), nlpar = "a"),
  prior(normal(0, 10), nlpar = "s", lb = 0, ub = Inf)
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                chains = 4,
                iter = 20000,
                backend = 'cmdstanr')

prior_chains = get_chains(prior_fit, type = "logreg_one_nl")

(p_prior_ce = plot_ce(prior_fit, plot_data = NA, index = 2, title = "Prior Predictive"))
(p_priors = plot_chains(prior_chains, plot_data = NA, color = "orange", title = "Prior Distributions", show_pointinterval = F))


(p_prior_ce = plot_ce(prior_fit, NA, 2))
(p_priors = plot_chains(prior_chains, NA, "orange", "Prior Distributions", F))


# posterior fit
# ===============================================
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

posterior_chains = get_chains(posterior_fit, type = "logreg_one_nl")
pars = get_pars(posterior_chains, d_sim)

# tbl = pars %>%
#   ggtexttable(rows = NULL,
#               theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_posterior_ce = plot_ce(posterior_fit, plot_data = d_sim, index = 2, title = "Posterior Predictive"))
(p_posterior = plot_chains(posterior_chains, plot_data = d_sim, color = 'cyan', title = "Posterior Distributions", show_pointinterval = T))
# (p_combo = plot_prior_vs_posterior(prior_chains, posterior_chains))

# plot overall
# -----------------------------------------------
((p_priors / p_posterior) | (p_prior_ce / p_posterior_ce))

