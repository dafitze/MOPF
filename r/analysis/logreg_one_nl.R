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
                               seq(from = -2.14, to = 2.14, length.out = 50)
                               # c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100
)
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

####################
# REPEAT PROCEDURE
####################

initial_fit = brm(model,
                  prior = priors,
                  data = d_sim,
                  cores = parallel::detectCores(),
                  chains = 4,
                  iter = 2000,
                  backend = 'cmdstanr')


tmp = tibble(idx = 1:50,
             b0 = 0.0,
             b1_1 = 1,
             # b1_2 = 3,
             # b1_3 = 6,
             # b1_4 = 9
             ) %>%
  mutate(dat1 = map2(b0, b1_1, ~cell_mean_simulation(b0_pre = .x, b1_pre = .y, n_vpn = 1, n_trials = 20, time = c("pre"), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)),
         # dat2 = map2(b0, b1_2, ~cell_mean_simulation(b0_pre = .x, b1_pre = .y, n_vpn = 1, n_trials = 20, time = c("pre"), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)),
         # dat3 = map2(b0, b1_3, ~cell_mean_simulation(b0_pre = .x, b1_pre = .y, n_vpn = 1, n_trials = 20, time = c("pre"), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)),
         # dat4 = map2(b0, b1_4, ~cell_mean_simulation(b0_pre = .x, b1_pre = .y, n_vpn = 1, n_trials = 20, time = c("pre"), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)),
         
         fit1 = map(dat1, ~update(initial_fit, newdata = .x)),
         # fit2 = map(dat2, ~update(initial_fit, newdata = .x)),
         # fit3 = map(dat3, ~update(initial_fit, newdata = .x)),
         # fit4 = map(dat4, ~update(initial_fit, newdata = .x)),
         
         chains1 = map(fit1, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus)),
         # chains2 = map(fit2, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus)),
         # chains3 = map(fit3, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus)),
         # chains4 = map(fit4, ~spread_draws(.x, b_inter_Intercept, b_preds_stimulus))
  )

tmp2_sim = tibble(param = rep(c('Intercept', 'stimulus'), each = 4),
                 nb = rep(c('1','2','3','4'), times = 2),
                 value = c(rep(unique(tmp$b0), 4), unique(tmp$b1_1), unique(tmp$b1_2), unique(tmp$b1_3), unique(tmp$b1_4)))
tmp_sim = tibble(param = c('b_inter_Intercept', 'b_preds_stimulus'),
                 value = c(unique(tmp$b0), unique(tmp$b1_1)))
tmp2 = bind_cols(
  chain1 = select(bind_rows(tmp$chains1, .id = "id"), id, b_inter_Intercept, b_preds_stimulus),
  chain2 = select(bind_rows(tmp$chains2, .id = "id"), b_inter_Intercept, b_preds_stimulus),
  chain3 = select(bind_rows(tmp$chains3, .id = "id"), b_inter_Intercept, b_preds_stimulus),
  chain4 = select(bind_rows(tmp$chains4, .id = "id"), b_inter_Intercept, b_preds_stimulus)
) |>
  pivot_longer(cols = -id, names_to = "param", values_to = "value") |>
  separate(param, into = c('class', 'par', 'param', 'nb'))

tmp$chains1 |>
  bind_rows(.id = "id") |>
  select(id, b_inter_Intercept, b_preds_stimulus) |>
  pivot_longer(cols = -id, names_to = "param", values_to = "value") |>
  ggplot(aes(x = value, y = id, color = id, xintercept = value)) +
  stat_halfeye(.width = c(.51,.93), geom = "pointinterval", alpha = 0.6, show.legend = F) +
  geom_vline(data = tmp_sim, aes(xintercept = value)) +
  facet_grid(~param, scales = "free") +
  labs(x = '',
       y = '') +
  theme_clean() +
  theme(
    # axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 

         