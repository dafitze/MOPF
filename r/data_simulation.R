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