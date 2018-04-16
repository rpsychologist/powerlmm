## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(powerlmm)

## ------------------------------------------------------------------------
p1 <- study_parameters(n1 = 11,
                      n2 = 25,
                      sigma_subject_intercept = 1.44,
                      sigma_subject_slope = 0.2,
                      sigma_error = 1.44,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))

p2 <- study_parameters(n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))
p1

## ------------------------------------------------------------------------
get_power(p2)

## ------------------------------------------------------------------------
p2 <- study_parameters(n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      dropout = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))

## ---- fig.width=8--------------------------------------------------------
plot(p2)

## ---- message = FALSE----------------------------------------------------
get_power(p2)

## ---- fig.width=4--------------------------------------------------------

d <- per_treatment(control = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                   treatment = dropout_weibull(proportion = 0.5, 
                                                rate = 2))

p2 <- study_parameters(n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      dropout = d,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))

plot(p2, type = "dropout")

## ------------------------------------------------------------------------
p1 <- study_parameters(n1 = 11,
                      n2 = 30,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))
p2 <- study_parameters(n1 = 11,
                      n2 = per_treatment(control = 10,
                                         treatment = 50),
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))

p1
p2

## ------------------------------------------------------------------------
get_power(p1)$power
get_power(p2)$power

## ------------------------------------------------------------------------
p1 <- study_parameters(n1 = 11,
                      n2 = 30,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.019,
                      effect_size = cohend(-0.5, 
                                           standardizer = "pretest_SD"))

x <- get_power_table(p1, n2 = seq(10, 30, by = 5), var_ratio = c(0.01, 0.02, 0.05))
x

## ---- fig.width = 5------------------------------------------------------
plot(x)

