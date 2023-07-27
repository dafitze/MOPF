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
source("r/plot/plot_ce.R")
source("r/plot/plot_pf.R")
source("r/plot/plot_param_recov.R")
d_sim = simulate_data(b0 = 0.0,
                      b0_sigma = 0.8,
                      b1 = 2.0,
                      b1_sigma = 0.4,
                      b2 = 0.0,
                      b2_sigma = 0.2,
                      b3 = 2.0,
                      b3_sigma = 0.4,
                      rho = -0.7,
                      n_vpn = 6,
                      n_trials = 20,
                      time = c("pre","post"),
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_sim, mu = 0.0)

d_sim_summary = d_sim %>%
  group_by(vpn, time, stimulus) %>%
  summarise(mean_response = mean(response)) %>%
  mutate(lower__ = 0,
         upper__ = 0,
         effect2__ = time)


# model
# -----------------------------------------------
# model = bf(response ~ 0 + Intercept + stimulus + (0 + Intercept + stimulus | vpn),
           # family = bernoulli('logit'))

model = bf(
  response ~ 0 + Intercept + time:stimulus + (0 + Intercept + time:stimulus | vpn),
  family = bernoulli('logit')
)

# prior
# -----------------------------------------------
# get_prior(model, d_sim)

priors = c(
  prior(normal(0.0, 10), class = "b", coef = "Intercept"),
  prior(normal(0.0, 10), class = "b", coef = "timepost:stimulus"),
  prior(normal(0.0, 10), class = "b", coef = "timepre:stimulus"),  
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(3), class = "cor")
)

# priors |>
#   parse_dist(prior) |>
#   filter(class == "cor") |>
#   # filter(nlpar == "eta") |>
#   ggplot(aes(y = class, dist = .dist, args = .args)) +
#   stat_dist_halfeye() +
#   facet_wrap(~coef, scales = "free_y") +
#   labs(
#     title = "stat_dist_halfeyeh()",
#     subtitle = "with brms::prior() and tidybayes::parse_dist() to visualize priors",
#     x = NULL
#   ) +
#   theme_tidybayes()
# 
# data.frame(prior = "lkjcorr(4)") %>%
#   parse_dist(prior) %>%
#   marginalize_lkjcorr(K = 2) %>%
#   ggplot(aes(y = prior, xdist = .dist_obj)) +
#   stat_halfeye() +
#   xlim(-1, 1) +
#   xlab("Marginal correlation for LKJ(3) prior on 2x2 correlation matrix")

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

(p_posterior_ce = plot_ce(posterior_fit, NA, 1))

# parameter estimates 
# -----------------------------------------------
get_variables(posterior_fit)
chains = posterior_fit %>%
  spread_draws(
    b_Intercept,
    `b_timepre:stimulus`,
    `b_timepost:stimulus`,
    sd_vpn__Intercept,
    `sd_vpn__timepre:stimulus`,
    `sd_vpn__timepost:stimulus`,
    `cor_vpn__Intercept__timepre:stimulus`,
    `cor_vpn__Intercept__timepost:stimulus`,
    `cor_vpn__timepre:stimulus__timepost:stimulus`
  ) %>%
  mutate(
    b0 = b_Intercept,
    b1_pre = `b_timepre:stimulus`,
    b1_post = `b_timepost:stimulus`,
    sd_b0 = sd_vpn__Intercept,
    sd_b1_pre = `sd_vpn__timepre:stimulus`,
    sd_b1_post = `sd_vpn__timepost:stimulus`,
    cor_b0_b1_pre = `cor_vpn__Intercept__timepre:stimulus`,
    cor_b0_b1_post = `cor_vpn__Intercept__timepost:stimulus`,
    cor_b1_prepost = `cor_vpn__timepre:stimulus__timepost:stimulus`,
    b1_diff = b1_post - b1_pre
    # mu = -b0 / b1,
    # sigma = 1 / b1
  ) %>%
  select(b0, b1_pre, b1_post, sd_b0, sd_b1_pre, sd_b1_post, cor_b0_b1_pre, cor_b0_b1_post, cor_b1_prepost, b1_diff)

pars = chains %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>% 
  group_by(param) %>%
  mean_qi(.width = 0.93) %>%
  select(param, value, .lower, .upper, .width) %>%
  arrange(factor(param, levels = c("b0", "sd_b0", "b1_pre", "sd_b1_pre", "b1_post", "sd_b1_post", "b1_diff", "cor_b0_b1_pre", "cor_b0_b1_post", "cor_b1_prepost"))) %>%
  mutate(fit = round(value, digits = 2),
         .lower = round(.lower, digits = 2),
         .upper = round(.upper, digits = 2),
         sim = c(
           unique(d_sim$b0),
           unique(d_sim$b0_sigma),
           unique(d_sim$b1),
           unique(d_sim$b1_sigma),
           unique(d_sim$b1) + unique(d_sim$b3),
           unique(d_sim$b1_sigma),
           unique(d_sim$b3),
           unique(d_sim$rho),
           unique(d_sim$rho),
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
  