p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.03)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f3 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (1 + time | cluster)")
f <- sim_formula_compare("m0" = f0, "m1" = f1, "m2" = f2, "m3" = f3)

res <- simulate(p, formula = f, nsim = 500, satterthwaite = TRUE, cores = 10, CI = FALSE)

# type 1 error increased
summary(res, model_selection = "FW")

# more liberal selection,
# type 1 error now 0.07
summary(res, model_selection = "FW", LRT_alpha = 0.25)

# compare with the correct model
summary(res, model = "m2")

# unecessary 3-level random slope
# conservative, and convergence warnings.
# leads overestiamed cluster-level random slope
summary(res, model = "m3")


res1 <- simulate(p, formula = f0, nsim = 100, satterthwaite = TRUE, cores = 10, CI = FALSE)



## multi
p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = c(0.01, 0.03))

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)
res_m <- simulate(p, formula = f, nsim = 1000, satterthwaite = TRUE, cores = 10, CI = FALSE)

x <- summary(res_m, model = "m2")
