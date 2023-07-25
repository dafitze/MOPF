library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(ggpmisc)
library(ggpubr)

# ===============================================
# Simulated Data
# ===============================================
source("r/data_simulation.R")
source("r/plot_ce.R")
source("r/param_recov_plot.R")
d_sim = simulate_ranef(
  b0 = 0.0,
  b1 = 2.0, 
  guess = 0,
  lapse = 0,
  sigma_b0 = 1,
  sigma_b1 = 0.8,
  rho = -0.7,
  n_vpn = 10,
  n_trials = 20
)


# plot_pf(d_sim, mu = 0.5)

d_sim_summary = d_sim %>%
  group_by(stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0)


# model
# -----------------------------------------------
# model = bf(response ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
           # family = bernoulli('logit'))

model = bf(
  response ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
  family = bernoulli('logit')
)

# prior
# -----------------------------------------------
# get_prior(model, d_sim)

priors = c(
  prior(normal(0.0, 10), class = "b", coef = "Intercept"),
  prior(normal(0.0, 10), class = "b", coef = "stimulus"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(3), class = "cor")
)

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = 'cmdstanr')

(p_prior_ce = plot_ce(prior_fit, NA, 1))


# posterior fit
# -----------------------------------------------
posterior_fit = brm(model,
                    prior = priors,
                    data = d_sim,
                    cores = parallel::detectCores(),
                    chains = 4,
                    iter = 2000,
                    backend = 'cmdstanr')

(p_posterior_ce = plot_ce(posterior_fit, d_sim_summary, 1))


# parameter estimates
# -----------------------------------------------
get_variables(posterior_fit)
chains = posterior_fit %>%
  spread_draws(
    b_Intercept,
    b_stimulus,
    sd_vpn__Intercept,
    sd_vpn__stimulus,
    cor_vpn__Intercept__stimulus
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1 = b_stimulus,
    sd_b0 = sd_vpn__Intercept,
    sd_b1 = sd_vpn__stimulus,
    cor_b0_b1 = cor_vpn__Intercept__stimulus
  ) %>%
  select(b0, b1, sd_b0, sd_b1, cor_b0_b1)


  pars = chains %>%
    pivot_longer(cols = 1:5, names_to = "param", values_to = "value") %>% 
    group_by(param) %>%
    mean_qi() %>%
    select(param, value, .lower, .upper, .width) %>%
    mutate(fit = round(value, digits = 2),
           .lower = round(.lower, digits = 2),
           .upper = round(.upper, digits = 2),
           sim = c(
             unique(d_sim$b0),
             unique(d_sim$b1),
             unique(d_sim$rho),
             unique(d_sim$sigma_b0),
             unique(d_sim$sigma_b1)
         )) %>%
  select(param, sim, fit, .lower, .upper)

tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
p_param_recov = ggplot(data = pars) +
  geom_pointinterval(mapping = aes(y = 1, x = fit, xmin = .lower, xmax = .upper, color = "95% CI")) +
  geom_point(mapping = aes(x = sim, y = 1, color = "Simulation"), 
             shape = 124, 
             size = 8) +
  facet_wrap(~param, 
             scales = "free_x",
             ncol = 1) +
  labs(x = "",
       y = "",
       color = "alsdkfj") +
  scale_color_manual(breaks = c("Simulation", "95% CI"),
                     values = c('Simulation'='red', '95% CI'='black')) +
  guides(color = guide_legend(title = '',
                              override.aes=list(shape = ""))) +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom")


(p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)


# individual values
# -----------------------------------------------
Intercept = coef(posterior_fit, pars = "Intercept")[[1]] %>%
  as_tibble() %>%
  select(1)
slope = coef(posterior_fit, pars = "stimulus")[[1]] %>%
  as_tibble() %>%
  select(1)
vp_pars = cbind(Intercept, slope)

tibble(Intercept = coef(posterior_fit, pars = "Intercept")[[1]]$estimate,
       b1 = coef(posterior_fit, pars = "stimulus")[[1]])
  
