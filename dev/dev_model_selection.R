p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0,
                      var_ratio = 0.03)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)


res <- simulate(p, formula = f, nsim = 100, satterthwaite = TRUE, cores = 10, CI = FALSE)
summary(res)

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
