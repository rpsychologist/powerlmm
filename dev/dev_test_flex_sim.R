

## simulate

f <- sim_formula("y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)")
f <- sim_formula("y ~ time + treatment + time:treatment + (1 + time | subject) + (0 + treatment:time | cluster)", test = c("time:treatment"))
f <- sim_formula("y ~ time + treatment + treatment:time + (1 + time | subject) + (0 + time | cluster)")
f1 <- sim_formula("y ~ treatment + (0 + treatment | cluster)", data_transform = transform_to_posttest, test = "treatment")

res <- simulate(p, formula = f, nsim = 50, satterthwaite = TRUE, cores = 1, CI = FALSE)
summary(res)

res <- simulate(p, nsim = 50, satterthwaite = TRUE)
res1 <- simulate(p, formula = f1, nsim = 10, CI = TRUE)

res1
summary(res1)

## compare
f2 <- compare_sim_formulas("3-lvl" = f, "2-lvl" = f1)
res2 <- simulate(p, formula = f2, nsim = 2000, cores = 15, satterthwaite = TRUE)
summary(res2)

## TODO: summary(res, para = "treatment")

p <- update(p, effect_size = cohend(1))


f <- sim_formula("y ~ treatment + pretest + (1 | cluster)", data_transform = transform_to_posttest, satterth_test = "treatment")
f2 <- compare_sim_formulas("3-lvl" = f, "2-lvl" = f1)
res2 <- simulate(p, formula = f2, nsim = 500, cores = 10, satterthwaite = TRUE)
summary(res2)

## TODO: fix so summary show all terms, e.g. 'pretest'




res3 <- simulate(update(p, n2 = c(5,10)), formula = f2, nsim = 5, cores = 1, satterthwaite = TRUE, CI = TRUE)
res3

summary(res3, para = "treatment")


f2 <- compare_sim_formulas("sim1" = f, "sim2" = f1, "sim3" = f1)
res4 <- simulate(update(p, n2 = c(5,10)), formula = f2, nsim = 5, cores = 1, satterthwaite = TRUE)

res4


### TODO: test with lmerTest < 3.0.0


## Stepwise

f <- sim_formula("y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)")
res <- simulate(p, nsim = 10, formula = f, satterthwaite = TRUE)
summary(res)


f0 <- sim_formula("y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)")
f1 <- sim_formula("y ~ treatment + (1 + time | subject)")
f <- compare_sim_formulas("3lvl" = f0, "2lvl" = f1)


res <- simulate(p, formula = f, nsim = 5, satterthwaite = TRUE, cores = 1, CI = FALSE)
summary(res)

res <- simulate(p, nsim = 50, satterthwaite = TRUE)



###


p <- study_parameters(n1 = 11,
                      n2 = 5,
                      n3 = 4,
                      fixed_slope = -2,
                      fixed_intercept = 1.2,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      dropout = dropout_weibull(0.3, 0.3),
                      effect_size = -1)


f <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)", test = "time1")
f <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)", test = c("time","treatment:time"))
f <- sim_formula("y ~ treatment * time + I(time^2) + (1 + time | subject) + (0 + time | cluster)")
res <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE, CI = TRUE,
                progress = FALSE)

res <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE, CI = TRUE,
                progress = FALSE)
summary(res)


## multi

p <- study_parameters(n1 = 11,
                      n2 = 5:6,
                      n3 = 4,
                      fixed_slope = -2,
                      fixed_intercept = 1.2,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      dropout = dropout_weibull(0.3, 0.3),
                      effect_size = -1)

f <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)", test = "treatment:time")
f <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)", test = "time:treatment")
res0 <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE, CI = TRUE,
                 progress = FALSE)
summary(res0)


## step
f0 <- sim_formula("y ~ time * treatment + (1 | subject)", test = "treatment:time")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)", test = "treatment:time")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)", test = "treatment:time")

res <- simulate(p, nsim = 2, formula = compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2), satterthwaite = TRUE, CI = TRUE,
                progress = FALSE)
summary(res, model_selection = "FW")

## transform

f0 <- sim_formula("y ~ treatment + (1 | cluster)", test = "treatment", data_transform = transform_to_posttest)

res <- simulate(p, nsim = 500, formula = f0, satterthwaite = TRUE, CI = FALSE, cores = 15,
                progress = FALSE)
summary(res, para = "treatment")


f0 <- sim_formula("y ~ 1 + (1 | cluster)", test = "treatment", data_transform = transform_to_posttest)
f1 <- sim_formula("y ~ treatment + (1 | cluster)", test = "treatment", data_transform = transform_to_posttest)

compare_sim_formulas("m0" = f0, "f1" = f1)



