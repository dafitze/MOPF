simulate_data = function(b0 = 0.0,
                         b0_sigma = 0.0,
                         b1 = 0.0,
                         b1_sigma = 0.0,
                         b2 = 0.0,
                         b2_sigma = 0.0,
                         b3 = 0.0,
                         b3_sigma = 0.0,
                         b0_guess = 0.0,
                         b0_guess_sigma = 0.0,
                         b0_lapse = 0.0,
                         b0_lapse_sigma = 0.0,
                         rho = 0.0,
                         n_vpn = 2,
                         n_trials = 20,
                         time = c("pre","post"),
                         stimulus = c(-2, -1,1,2)){
  
  
  if (n_vpn == 1){                                                  # SIMULATE DATA FOR ONE SUBJECT
    vpn_param = tibble(vpn = 1, b0, b1, b0_lapse, b0_guess)
  } else {                                                          # SIMULATE DATA FOR MULTIPLE SUBJECTS
    mus = c(b0, b1, b0_lapse, b0_guess)                             #   means
    sigmas = c(b0_sigma, b1_sigma, b0_lapse_sigma, b0_guess_sigma)  #   standard deviations
    rho_ma = matrix(c(1, rho, rho, rho,                             #   correlation matrix
                      rho, 1, rho, rho,
                      rho, rho, 1, rho,
                      rho, rho, rho, 1), nrow = 4)
    SD = diag(sigmas) %*% rho_ma %*% diag(sigmas)                   #   now matrix multiply to get covariance matrix
    
    vpn_param = MASS::mvrnorm(n_vpn, mus, SD, tol = 0.1) |>         #   draw param
      as_tibble(.name_repair = "unique") |>
      purrr::set_names("b0", "b1", "b0_lapse", "b0_guess") |>
      mutate(vpn = 1:n_vpn, .before = "b0")
  }  
  
  df = vpn_param |>                                                # THE SAME FOR ONE OR MULTIPLE
    expand_grid(time, 
                stimulus,
                rep = seq(1, n_trials)) |> 
    mutate(
      time = factor(time, levels = c("pre", "post")),
      time_num = case_when(time == "pre" ~ 0,
                           time == "post" ~ 1),
      guess = b0_guess,# + (b2_guess * time_num),
      lapse = b0_lapse,# + (b2_lapse * time_num),
      linpred = b0 + (b1 * stimulus) + (b2 * time_num) + (b3 * stimulus * time_num),
      theta = guess + (1 - guess - lapse) * plogis(linpred),
      response = rbinom(n = length(theta), size = 1, prob = theta),
    ) 
  # |>
    # arrange(vpn, stimulus)
}






# OLD FUNCTIONS
# --------------
# 
# simulate_data = function(b0 = 0.0,
#                          b1 = 2.0, 
#                          guess = 0,
#                          lapse = 0,
#                          vpn = 1,
#                          nreps = 20) {
#   
#   df = expand_grid(
#     vpn = vpn,
#     rep = seq(1, nreps),
#     stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
#     b0 = b0,
#     b1 = b1,
#     guess = guess,
#     lapse = lapse
#   ) %>%
#     mutate(theta = guess + (1 - guess - lapse) * plogis(b0 + b1 * stimulus),
#            response = rbinom(n = length(theta), size = 1, prob = theta)
#     )
#   df <- df %>%
#     arrange(vpn, stimulus, rep)
#   return(df)
# }
# 
# simulate_prepost_data = function(b0 = 0.0,
#                                  b1 = 2.0, 
#                                  b2 = 0.1,
#                                  b3 = 3.0,
#                                  b0_guess = 0.0,
#                                  b2_guess = 0.0,
#                                  b0_lapse = 0.0,
#                                  b2_lapse = 0.0,
#                                  vpn = 1,
#                                  nreps = 20,
#                                  time = c("pre","post")) {
#   
#   df = expand_grid(
#     vpn = vpn,
#     time = time,
#     rep = seq(1, nreps),
#     stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100,
#     b0 = b0,
#     b1 = b1,
#     b2 = b2,
#     b3 = b3,
#     b1_pre = b1,
#     b1_post = b1 + b3,
#     b0_guess = b0_guess,
#     b2_guess = b2_guess,
#     b0_lapse = b0_lapse,
#     b2_lapse = b2_lapse
#   ) %>%
#     mutate(
#       time = factor(time, levels = c("pre", "post")),
#       time_num = case_when(time == "pre" ~ 0,
#                            time == "post" ~ 1),
#       guess = b0_guess + (b2_guess * time_num),
#       lapse = b0_lapse + (b2_lapse * time_num),
#       theta = guess + (1 - guess - lapse) * plogis(b0 + (b1 * stimulus) + (b2 * time_num) + (b3 * stimulus * time_num)),
#       response = rbinom(n = length(theta), size = 1, prob = theta)
#     )
#   df <- df %>%
#     arrange(vpn, stimulus, rep)
#   return(df)
# }
# 
# simulate_one_ranef = function(b0 = 0.0,
#                               b1 = 2.0, 
#                               guess = 0,
#                               lapse = 0,
#                               sigma_b0 = 1,
#                               sigma_b1 = 0.8,
#                               rho = -0.7,
#                               n_vpn = 5,
#                               n_trials = 20) {
#   
#   # Ccombine the terms
#   mu <- c(b0, b1)
#   cov_b01 <- sigma_b0 * sigma_b1 * rho
#   SD  <- matrix(c(sigma_b0^2, cov_b01,
#                   cov_b01, sigma_b1^2), ncol = 2)
#   
#   varying_effects <-
#     MASS::mvrnorm(n_vpn, mu, SD) |>
#     # as_tibble(.name_repair = "unique") |>
#     data.frame() |>
#     purrr::set_names("b0_j", "b1_j")
#   
#   d <- varying_effects |>
#     mutate(vpn  = 1:n_vpn) |>
#     expand(nesting(vpn, b0_j, b1_j), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |>
#     mutate(linpred = b0_j + b1_j * stimulus) |>
#     slice(rep(1:n(), each = n_trials)) |>
#     mutate(theta = plogis(linpred),
#            response = rbinom(n = length(theta), size = 1, prob = theta),
#            b0 = b0,
#            b1 = b1,
#            sigma_b0 = sigma_b0,
#            sigma_b1 = sigma_b1,
#            rho = rho,)
#   
#   d
# }
# 
# simulate_prepost_ranef = function(b0 = 0.0,
#                                   b1 = 2.0,
#                                   b2 = 0.1,
#                                   b3 = 0.3,
#                                   guess = 0,
#                                   lapse = 0,
#                                   sigma_b0 = 1,
#                                   sigma_b1 = 0.8,
#                                   sigma_b2 = 0.8,
#                                   sigma_b3 = 0.8,
#                                   rho = -0.7,
#                                   n_vpn = 3,
#                                   n_trials = 2,
#                                   sigma = 0.5){
#   
#   # Ccombine the terms
#   mu <- c(b0, b1, b2, b3)
#   cov_b01 = sigma_b0 * sigma_b1 * rho 
#   cov_b02 = sigma_b0 * sigma_b2 * rho
#   cov_b03 = sigma_b0 * sigma_b3 * rho
#   cov_b12 = sigma_b1 * sigma_b2 * rho
#   cov_b13 = sigma_b1 * sigma_b3 * rho
#   cov_b23 = sigma_b2 * sigma_b3 * rho
#   SD  <- matrix(c(sigma_b0^2, cov_b01, cov_b02, cov_b03,
#                   cov_b01, sigma_b1^2, cov_b12, cov_b13,
#                   cov_b02, cov_b12, sigma_b2^2, cov_b23,
#                   cov_b03, cov_b13, cov_b23, sigma_b3^2), 
#                 ncol = 4)
#   
#   #  sigmas <- c(sigma_a, sigma_b)          # standard deviations
#   #  rho <- matrix(c(1, rho,             # correlation matrix
#   #                 rho, 1), nrow = 2)
#   #
#   # # now matrix multiply to get covariance matrix
#   # SD <- diag(sigmas) %*% rho %*% diag(sigmas)
#   
#   varying_effects <-
#     MASS::mvrnorm(n_vpn, mu, SD, tol = 1) |>
#     # as_tibble(.name_repair = "unique") |>
#     data.frame() |>
#     purrr::set_names("b0_j", "b1_j", "b2_j", "b3_j")
#   
#   d <- varying_effects |>
#     mutate(vpn  = 1:n_vpn) |>
#     expand(nesting(vpn, b0_j, b1_j, b2_j, b3_j), time_num = c(0, 1), stimulus = c(-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214)/100) |>
#     mutate(linpred = b0_j + (b1_j * stimulus) + (b2_j * time_num) + (b3_j * stimulus * time_num)) |>
#     mutate(time = ifelse(time_num == 0, "pre", "post"),
#            time = factor(time, levels = c("pre", "post"))) |>
#     slice(rep(1:n(), each = n_trials)) |>
#     mutate(theta = plogis(linpred),
#            response = rbinom(n = length(theta), size = 1, prob = theta),
#            b0 = b0,
#            b1 = b1,
#            b2 = b2,
#            b3 = b3,
#            b1_pre = b1,
#            b1_post = b1 + b3,
#            sigma_b0 = sigma_b0,
#            sigma_b1 = sigma_b1,
#            sigma_b2 = sigma_b2,
#            sigma_b3 = sigma_b3,
#            rho = rho)
# 
#   d
# }
# 




