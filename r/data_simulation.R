simulate_data = function(b0 = 0.0,
                         b1 = 2.0, 
                         guess = 0,
                         lapse = 0,
                         vpn = 1,
                         nreps = 20) {
  
  df = expand_grid(
    vpn = vpn,
    rep = seq(1, nreps),
    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
    b0 = b0,
    b1 = b1,
    guess = guess,
    lapse = lapse
  ) %>%
    mutate(theta = guess + (1 - guess - lapse) * plogis(b0 + b1 * stimulus),
           response = rbinom(n = length(theta), size = 1, prob = theta)
    )
  df <- df %>%
    arrange(vpn, stimulus, rep)
  return(df)
}

simulate_prepost_data = function(b0 = 0.0,
                                 b1 = 2.0, 
                                 b2 = 0.1,
                                 b3 = 0.3,
                                 guess = 0,
                                 lapse = 0,
                                 vpn = 1,
                                 nreps = 20,
                                 time = c("pre","post")) {
  
  df = expand_grid(
    vpn = vpn,
    time = time,
    rep = seq(1, nreps),
    stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b1_pre = b1,
    b1_post = b1 + b3,
    guess = guess,
    lapse = lapse
  ) %>%
    mutate(
      time = factor(time, levels = c("pre", "post")),
      time_num = case_when(time == "pre" ~ 0,
                           time == "post" ~ 1),
      theta = guess + (1 - guess - lapse) * plogis(b0 + (b1 * stimulus) + (b2 * time_num) + (b3 * stimulus * time_num)),
      response = rbinom(n = length(theta), size = 1, prob = theta)
    )
  df <- df %>%
    arrange(vpn, stimulus, rep)
  return(df)
}

simulate_one_ranef = function(b0 = 0.0,
                              b1 = 2.0, 
                              guess = 0,
                              lapse = 0,
                              sigma_b0 = 1,
                              sigma_b1 = 0.8,
                              rho = -0.7,
                              n_vpn = 5,
                              n_trials = 20) {
  
  # Ccombine the terms
  mu <- c(b0, b1)
  cov_b01 <- sigma_b0 * sigma_b1 * rho
  SD  <- matrix(c(sigma_b0^2, cov_b01,
                  cov_b01, sigma_b1^2), ncol = 2)
  
  varying_effects <-
    MASS::mvrnorm(n_vpn, mu, SD) |>
    # as_tibble(.name_repair = "unique") |>
    data.frame() |>
    purrr::set_names("b0_j", "b1_j")
  
  d <- varying_effects |>
    mutate(vpn  = 1:n_vpn) |>
    expand(nesting(vpn, b0_j, b1_j), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |>
    mutate(linpred = b0_j + b1_j * stimulus) |>
    slice(rep(1:n(), each = n_trials)) |>
    mutate(theta = plogis(linpred),
           response = rbinom(n = length(theta), size = 1, prob = theta),
           b0 = 0.0,
           b1 = 2.0,
           sigma_b0 = 1,
           sigma_b1 = 0.8,
           rho = -0.7,)
  
  d
}

simulate_prepost_ranef = function(b0 = 0.0,
                                  b1 = 2.0,
                                  b2 = 0.1,
                                  b3 = 0.3,
                                  guess = 0,
                                  lapse = 0,
                                  sigma_b0 = 1,
                                  sigma_b1 = 0.8,
                                  sigma_b2 = 0.8,
                                  sigma_b3 = 0.8,
                                  rho = -0.7,
                                  n_subjects = 3,
                                  n_trials = 2,
                                  sigma = 0.5){
  
  # Ccombine the terms
  mu <- c(b0, b1, b2, b3)
  cov_b01 = sigma_b0 * sigma_b1 * rho 
  cov_b02 = sigma_b0 * sigma_b2 * rho
  cov_b03 = sigma_b0 * sigma_b3 * rho
  cov_b12 = sigma_b1 * sigma_b2 * rho
  cov_b13 = sigma_b1 * sigma_b3 * rho
  cov_b23 = sigma_b2 * sigma_b3 * rho
  SD  <- matrix(c(sigma_b0^2, cov_b01, cov_b02, cov_b03,
                  cov_b01, sigma_b1^2, cov_b12, cov_b13,
                  cov_b02, cov_b12, sigma_b2^2, cov_b23,
                  cov_b03, cov_b13, cov_b23, sigma_b3^2), 
                ncol = 4)
  
  #  sigmas <- c(sigma_a, sigma_b)          # standard deviations
  #  rho <- matrix(c(1, rho,             # correlation matrix
  #                 rho, 1), nrow = 2)
  #
  # # now matrix multiply to get covariance matrix
  # SD <- diag(sigmas) %*% rho %*% diag(sigmas)
  
  varying_effects <-
    MASS::mvrnorm(n_subjects, mu, SD, tol = 1) |>
    # as_tibble(.name_repair = "unique") |>
    data.frame() |>
    purrr::set_names("b0_j", "b1_j", "b2_j", "b3_j")
  
  d <- varying_effects |>
    mutate(vpn  = 1:n_subjects) |>
    expand(nesting(vpn, b0_j, b1_j), time = c(0, 1), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |>
    mutate(linpred = b0_j + b1_j * time) |>
    mutate(time_num = ifelse(time == 0, "pre", "post"),
           time_num = factor(time_num, levels = c("pre", "post"))) |>
    slice(rep(1:n(), each = n_trials)) |>
    mutate(theta = plogis(b0 + b1 * stimulus),
           response = rbinom(n = length(theta), size = 1, prob = theta))
  
  d
}  


