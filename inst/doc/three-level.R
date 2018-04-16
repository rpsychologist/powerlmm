## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(powerlmm)

## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
p

## ------------------------------------------------------------------------
get_power(p)

## ------------------------------------------------------------------------
get_VPC(p)

## ------------------------------------------------------------------------
get_correlation_matrix(p)

## ------------------------------------------------------------------------
get_sds(p)

## ---- fig.width = 8, message = FALSE, warning = FALSE--------------------
library(ggplot2)
x <- get_power_table(p, n2 = 5:20, n3 = c(4, 6, 8, 12), icc_slope = c(0.01, 0.05, 0.1))
plot(x) + scale_x_continuous(breaks = seq(20, 240, length.out = 5))

## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = per_treatment(control = 2, 
                                         treatment = 10),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD")
                      )
p

## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = per_treatment(control = 10,
                                         treatment = 2),
                      n3 = per_treatment(control = 2, 
                                         treatment = 10),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
p

## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(2, 5, 10, 30),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
p

## ------------------------------------------------------------------------

n2 <- per_treatment(control = unequal_clusters(5, 10, 15),
                    treatment = unequal_clusters(2, 3, 5, 5, 10, 15, 25))

p <- study_parameters(n1 = 11,
                      n2 = n2,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
p

## ------------------------------------------------------------------------
n2 <- unequal_clusters(func = rpois(n = 5, lambda = 5))

p <- study_parameters(n1 = 3,
                      n2 = n2,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
get_power(p, R = 10, progress = FALSE)

## ------------------------------------------------------------------------
# sample cluster sizes in each treatment group independently
# but from the same distribution
func <- unequal_clusters(func = rpois(n = 5, lambda = 5))
n2 <- per_treatment(control = func, 
                    treatment = func)


## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(2, 5, 10, 30),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      partially_nested = TRUE,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

p

## ------------------------------------------------------------------------
p1 <- study_parameters(n1 = 11,
                      n2 = 5,
                      n3 = 5,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      partially_nested = TRUE,                      
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

p2 <- study_parameters(n1 = 11,
                      n2 = per_treatment(control = 50, 
                                         treatment = 5),
                      n3 = per_treatment(control = 1,
                                         treatment = 5),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      partially_nested = TRUE,                      
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

p1
p2

## ------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 5,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0, 
                      var_ratio = 0.02,
                      icc_slope = 0.05,
                      dropout = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                      fixed_slope = -0.5/10,                      
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

p

## ---- fig.width=8--------------------------------------------------------
plot(p)

## ---- message = FALSE----------------------------------------------------
get_power(p)

## ---- fig.width=4--------------------------------------------------------
d <- per_treatment(control = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                   treatment = dropout_weibull(proportion = 0.5, 
                                                rate = 2))

p2 <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 5,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0, 
                      var_ratio = 0.02,
                      icc_slope = 0.05,
                      dropout = d,
                      fixed_slope = -0.5/10,                      
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

plot(p2, type = "dropout")

## ---- message = FALSE----------------------------------------------------
p2 <- study_parameters(n1 = 11,
                      n2 = c(5, 10, 15, 20, 30),
                      n3 = 5,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0, 
                      var_ratio = 0.02,
                      icc_slope = 0.05,
                      dropout = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                      fixed_slope = -0.5/10,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
get_DEFT(p2)

