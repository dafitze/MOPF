library(tidyverse)
library(manipulate)
library(ggplot2)
library(truncnorm)
library(crch) # truncated student t

# Psychometric Function
manipulate({
  ggplot(tibble(x = c(-2.5,2.5)), 
         aes(x)) +
    stat_function(fun = function(x) plogis(b0 + b1 * x), n = 1000) +
    geom_point(aes(x = -b0/b1, y = 0.5),
               shape = 124,
               color = 'red') +
    geom_vline(xintercept = -b0/b1, linetype = 2, alpha = .3) +
    geom_hline(yintercept = 0.5, linetype = 2, alpha = .3) +
    annotate("text", x = -1.5, y = 0.75, label = paste("mu = ", round(-b0/b1, digits = 3))) +
    annotate("text", x = -0.5, y = 0.75, label = paste("b0 = ", b0)) +
    annotate("text", x = -1.5, y = 0.7, label = paste("sigma = ", round(1/b1, digits = 3))) +
    annotate("text", x = -0.5, y = 0.7, label = paste("b1 = ", b1)) +
    lims(y = c(0,1)) +
    theme_minimal()
},
b0 = slider(min = -5, max = 5, initial = 0, step = 0.1),
b1 = slider(min = 0, max = 50, initial = 1)
)

# Prior Intercept
manipulate({
  ggplot(tibble(x = c(-25,25)), 
         aes(x)) +
    stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = mysd)) +
    geom_segment(aes(x = qnorm(c(0.025,0.975), mean = 0, sd = mysd), xend = qnorm(c(0.025,0.975), mean = 0, sd = mysd), 
                     y = 0, yend = Inf), linetype = 2, alpha = .3) +
    annotate("text", x = qnorm(0.025, mean = 0, sd = mysd) - 5, y = -0.01, label = paste("b0 = ", round(qnorm(0.025, mean = 0, sd = mysd),1))) +
    annotate("text", x = qnorm(0.975, mean = 0, sd = mysd) + 5, y = -0.01, label = paste("b0 = ", round(qnorm(0.975, mean = 0, sd = mysd),1))) +
    annotate("text", x = qnorm(0.025, mean = 0, sd = mysd) - 5, y = -0.02, label = paste("mu = ", round(-qnorm(0.025, mean = 0, sd = mysd)/ b1,1))) +
    annotate("text", x = qnorm(0.975, mean = 0, sd = mysd) + 5, y = -0.02, label = paste("mu = ", round(-qnorm(0.975, mean = 0, sd = mysd)/b1,1))) +
    ggtitle(paste("N(0,", mysd, ")", sep = "")) +
    theme_minimal()
},
mysd = slider(min = 0, max = 10, initial = 5),
b0 = slider(min = -5, max = 5, initial = 0, step = 0.1),
b1 = slider(min = 0, max = 50, initial = 1))

# Prior Stimulus
manipulate({
  ggplot(tibble(x = c(-45,45)), 
         aes(x)) +
    stat_function(fun = dtt, n = 1000, args = list(df = 1, location = 0, scale = mysd, left = 0, right = Inf)) +
    # geom_segment(aes(x = qtt(x, df = 1, location = 0, scale = mysd, df, left = 0, right = Inf, log = FALSE), xend = qnorm(0.975, mean = 0, sd = mysd),
                     # y = 0, yend = Inf), linetype = 2, alpha = .3) +
    annotate("text", x = qtruncnorm(0.975, a = 0, b = Inf, mean = 0, sd = mysd), y = -0.001, label = paste(round(qnorm(0.975, mean = 0, sd = mysd),1))) +
    ggtitle(paste("tt(0,", mysd, ")", sep = "")) +
    theme_minimal()
},
mysd = slider(min = 0, max = 20, initial = 10))

# Prior Lapse / Guess
manipulate({
  ggplot(tibble(x = c(0,1)), 
         aes(x)) +
    stat_function(fun = dbeta, n = 1000, args = list(shape1 = mya, shape2 = myb)) +
    # geom_segment(aes(x = qnorm(c(0.025,0.975), mean = 0, sd = mysd), xend = qnorm(c(0.025,0.975), mean = 0, sd = mysd), 
                     # y = 0, yend = Inf), linetype = 2, alpha = .3) +
    # annotate("text", x = qnorm(0.025, mean = 0, sd = mysd), y = -0.001, label = paste(round(qnorm(0.025, mean = 0, sd = mysd),1))) +
    # annotate("text", x = qnorm(0.975, mean = 0, sd = mysd), y = -0.001, label = paste(round(qnorm(0.975, mean = 0, sd = mysd),1))) +
    ggtitle(paste("Beta(", mya, ",", myb, ")", sep = "")) +
    theme_minimal()
},
mya = slider(min = 0, max = 50, initial = 2),
myb = slider(min = 0, max = 50, initial = 50)
)


# Prior Training

