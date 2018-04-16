## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(powerlmm)
knitr::opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  p <- study_parameters(n1 = 11,
#                        n2 = 10,
#                        n3 = 6,
#                        fixed_intercept = 37,
#                        fixed_slope = -0.64,
#                        sigma_subject_intercept = 2.8,
#                        sigma_subject_slope = 0.4,
#                        sigma_cluster_intercept = 0,
#                        cor_subject = -0.2,
#                        icc_slope = 0.05,
#                        sigma_error = 2.6,
#                        dropout = dropout_weibull(proportion = 0.3,
#                                                  rate = 1/2),
#                        partially_nested = TRUE,
#                        effect_size = cohend(-0.8,
#                                             standardizer = "pretest_SD"))
#  
#  p

## ------------------------------------------------------------------------
#  x <- get_power(p)
#  get_monte_carlo_se(x, nsim = c(1000, 2000, 5000))

## ------------------------------------------------------------------------
#  f <- "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)"
#  
#  cores <- parallel::detectCores() # use all cores
#  
#  res <- simulate(object = p,
#                  nsim = 5000,
#                  formula = f,
#                  satterthwaite = TRUE,
#                  cores = cores,
#                  save = FALSE)
#  
#  summary(res)

## ------------------------------------------------------------------------
#  p <- update(p, effect_size = 0)
#  
#  f <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)",
#            "wrong" = "y ~ treatment * time + (1 + time | subject)")
#  
#  res <- simulate(object = p,
#                  nsim = 5000,
#                  formula = f,
#                  satterthwaite = TRUE,
#                  cores = cores,
#                  save = FALSE)
#  
#  summary(res)

## ------------------------------------------------------------------------
#  p <- study_parameters(n1 = 11,
#                        n2 = 10,
#                        n3 = c(6, 12),
#                        fixed_intercept = 37,
#                        fixed_slope = -0.64,
#                        sigma_subject_intercept = 2.8,
#                        sigma_subject_slope = 0.4,
#                        sigma_cluster_intercept = 0,
#                        cor_subject = -0.2,
#                        icc_slope = c(0.05, 0.1),
#                        sigma_error = 2.6,
#                        dropout = dropout_weibull(proportion = 0.3,
#                                                  rate = 1/2),
#                        partially_nested = TRUE,
#                        effect_size = cohend(c(0.5, 0.8),
#                                             standardizer = "pretest_SD"))
#  
#  
#  f <- "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)"
#  
#  res <- simulate(object = p,
#                  nsim = 5000,
#                  formula = f,
#                  satterthwaite = TRUE,
#                  cores = cores,
#                  save = FALSE)

## ------------------------------------------------------------------------
#  # Summarize the 'time:treatment' results
#  summary(res, para = "time:treatment", type = "fixed", model = "correct")
#  
#  # Summarize the cluster-level random slope
#  summary(res, para = "cluster_slope", type = "random", model = "correct")
#  

