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
source("r/simulation/data_simulation.R")
source("r/plot/theme_clean.R")
source("r/plot/plot_ce.R")
source("r/plot/plot_pf.R")
source("r/plot/plot_priors.R")
source("r/plot/plot_param_recov.R")
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
  prior(normal(0, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 5), class = "b", coef = "stimulus"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(3), class = "cor")
)

(p_priors = plot_priors_logreg(priors))

prior_fit = brm(model,
                prior = priors,
                data = d_sim,
                sample_prior = "only",
                backend = 'cmdstanr')

(p_prior_ce = plot_ce(prior_fit, NA, 1))

prior_chains = prior_fit %>%
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
# get_variables(posterior_fit)
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
  mean_qi(.width = 0.93) %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(factor(param, levels = c("b0", "sd_b0", "b1", "sd_b1", "cor_b0_b1"))) %>%
  # arrange(match(param, c("b0", "sd_b0", "b1", "sd_b1", "cor_b0_b1"))) %>%
  mutate(fit = round(value, digits = 2),
         .lower = round(.lower, digits = 2),
         .upper = round(.upper, digits = 2),
         sim = c(
           unique(d_sim$b0),
           unique(d_sim$b0_sigma),
           unique(d_sim$b1),
           unique(d_sim$b1_sigma),
           unique(d_sim$rho)
         )) %>%
  select(param, sim, fit, .lower, .upper)

tbl = pars %>%
  # mutate(desciption = c("bias (logit-scale)", "slope (logit-scale)")) %>%
  ggtexttable(rows = NULL,
              theme = ttheme('blank'))

# parameter recovery
# -----------------------------------------------
(p_param_recov = plot_param_recov(chains, pars))


# plot overall
# -----------------------------------------------
(p_priors | p_prior_ce | p_posterior_ce) / (p_param_recov | tbl)





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
  
