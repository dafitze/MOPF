library(tidyverse)

source("r/data_simulation.R")
source("r/plot_pf.R")

# logreg one
# -----------------------------------------------
d_one = simulate_data(b0 = 0.0,
                      b1 = 2.0,
                      n_vpn = 1,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_one, mu = 0.0)

# logreg one (pre/post)
# -----------------------------------------------
d_one_prepost = simulate_data(b0 = 0.0,
                              b1 = 2.0,
                              b2 = 0.0,
                              b3 = 2.0,
                              n_vpn = 1,
                              n_trials = 20,
                              time = c("pre","post"),
                              stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_one_prepost, mu = 0.0)

# logreg multiple (pre/post)
# -----------------------------------------------
d_one_prepost = simulate_data(b0 = 0.0,
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

plot_pf(d_one_prepost, mu = 0.0)


# lapse one
# -----------------------------------------------
d_one = simulate_data(b0 = 0.0,
                      b1 = 5.0,
                      b0_guess = 0.1,
                      b0_lapse = 0.1,
                      n_vpn = 1,
                      n_trials = 20,
                      time = "pre",
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_one, mu = 0.0)

# lapse one (pre/post)
# -----------------------------------------------
d_one = simulate_data(b0 = 0.0,
                      b1 = 5.0,
                      b2 = 0.0,
                      b3 = 2.0,
                      b0_guess = 0.1,
                      b0_lapse = 0.1,
                      n_vpn = 1,
                      n_trials = 20,
                      time = c("pre","post"),
                      stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_one, mu = 0.0)


# lapse multiple (pre/post)
# -----------------------------------------------
d_one_prepost = simulate_data(b0 = 0.0,
                              b0_sigma = 0.8,
                              b1 = 2.0,
                              b1_sigma = 0.4,
                              b2 = 0.0,
                              b2_sigma = 0.2,
                              b3 = 2.0,
                              b3_sigma = 0.4,
                              b0_guess = 0.1,
                              b0_guess_sigma = 0.001,
                              b0_lapse = 0.1,
                              b0_lapse_sigma = 0.001,
                              rho = -0.7,
                              n_vpn = 6,
                              n_trials = 20,
                              time = c("pre","post"),
                              stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100)

plot_pf(d_one_prepost, mu = 0.0)
