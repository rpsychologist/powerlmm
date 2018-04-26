p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      partially_nested = TRUE,
                      var_ratio = 0.03)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)


res <- simulate(p, formula = f, nsim = 5, satterthwaite = TRUE, cores = 1, CI = FALSE)
summary(res)

winners <- step.plcp_sim(res)


# TODO: write function to use 'winners' to pick parameters from winning models and summarise
tmp <- lapply(winners, function(d) {
    ...
})


