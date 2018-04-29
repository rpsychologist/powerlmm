params <-
structure(list(EVAL = FALSE), .Names = "EVAL")

## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(powerlmm)
knitr::opts_chunk$set(
    eval = FALSE
    )

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
#  p

## ------------------------------------------------------------------------
#  x <- get_power(p)
#  get_monte_carlo_se(x, nsim = c(1000, 2000, 5000))

## ------------------------------------------------------------------------
#  f <- "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)"
#  
#  cores <- parallel::detectCores(logical = FALSE) # use all physical CPU cores
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
#  f <- sim_formula_compare("correct" = "y ~ treatment * time +
#                                        (1 + time | subject) +
#                                        (0 + treatment:time | cluster)",
#                           "wrong" = "y ~ treatment * time +
#                                      (1 + time | subject)")
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
#  f <- sim_formula("y ~ treatment * time + (1 + time | subject) +
#                   (0 + treatment:time | cluster)")
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

## ------------------------------------------------------------------------
#  cores <- parallel::detectCores(logical = FALSE)
#  
#  p <- study_parameters(n1 = 11,
#                        n2 = 40,
#                        n3 = 3,
#                        icc_pre_subject = 0.5,
#                        cor_subject = -0.5,
#                        icc_slope = 0.05,
#                        var_ratio = 0.03)
#  
#  f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
#  f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
#  f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
#  f3 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (1 + time | cluster)")
#  f <- sim_formula_compare("subject-intercept" = f0,
#                           "subject-slope" = f1,
#                           "cluster-slope" = f2,
#                           "cluster-intercept" = f3)
#  
#  res <- simulate(p, formula = f,
#                  nsim = 5000,
#                  satterthwaite = TRUE,
#                  cores = cores,
#                  CI = FALSE)
#  
#  # type 1 error increased
#  summary(res, model_selection = "FW")
#  
#  # but better then always ignoring
#  summary(res, model = "subject-slope")
#  
#  # more liberal selection,
#  # type 1 error now 0.07
#  summary(res, model_selection = "FW", LRT_alpha = 0.25)
#  
#  # compare with the correct model
#  summary(res, model = "m2")
#  
#  # unecessary 3-level random slope
#  # conservative, and convergence warnings.
#  # leads overestiamed cluster-level random slope
#  summary(res, model = "m3")
#  
#  
#  # See if power is impacted
#  p1 <- update(p, effect_size = cohend(0.8))
#  res_power <- simulate(p1,
#                        formula = f,
#                        nsim = 500,
#                        satterthwaite = TRUE,
#                        cores = cores,
#                        CI = FALSE)
#  
#  

## ------------------------------------------------------------------------
#  p <- study_parameters(n1 = 11,
#                        n2 = 20,
#                        n3 = 3,
#                        icc_pre_subject = 0.5,
#                        icc_pre_cluster = 0.1,
#                        icc_slope = 0.05,
#                        var_ratio = 0.03)
#  
#  # simulation formulas
#  # analyze as a posttest only 2-level model
#  f_pt <- sim_formula("y ~ treatment + (1 | cluster)",
#                   test = "treatment",
#                   data_transform = transform_to_posttest)
#  
#  # analyze as 3-level longitudinal
#  f_lt <- sim_formula("y ~ time*treatment +
#                           (1 + time | subject) +
#                           (1 + time | cluster)")
#  
#  f <- sim_formula_compare("posttest" = f_pt, "longitudinal" = f_lt)
#  ## Not run:
#  res <- simulate(p,
#                  formula = f,
#                  nsim = 2000,
#                  cores = cores,
#                  satterthwaite = TRUE)
#  summary(res)

## ------------------------------------------------------------------------
#  # simulation formulas
#  # analyze as a posttest only 2-level model
#  f0 <- sim_formula("y ~ treatment",
#                   test = "treatment",
#                   data_transform = transform_to_posttest)
#  
#  f1 <- sim_formula("y ~ treatment + (1 | cluster)",
#                   test = "treatment",
#                   data_transform = transform_to_posttest)
#  
#  f <- sim_formula_compare("post_ignore" = f0, "post_2lvl" = f1)
#  ## Not run:
#  res <- simulate(p,
#                  formula = f,
#                  nsim = 2000,
#                  cores = cores,
#                  satterthwaite = TRUE)
#  
#  # type I OLS model
#  sumamry(res, model = "post_ignore")
#  
#  # model selection
#  summary(res, model_selection = "FW")

