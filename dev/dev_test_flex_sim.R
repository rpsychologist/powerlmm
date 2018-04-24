f <- sim_formula("y ~ time + (1 | subject)")



transform_to_posttest(d)


p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.02)

d <- simulate_data(p)

##
check_formula(f)


f1 <- compare_sim_formulas("sim1" = f, "sim2" = sim_formula("y ~ time55 + (1 | subject)"))

str(f1)

names(f1)

check_formula(f1)


## simulate

f <- sim_formula("y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)")
f1 <- sim_formula("y ~ treatment + (1 | cluster)", data_transform = transform_to_posttest, satterth_test = "treatment")

res <- simulate(p, formula = f, nsim = 500, satterthwaite = TRUE, cores = 10)
res <- simulate(p, nsim = 50, satterthwaite = TRUE)
res1 <- simulate(p, formula = f1, nsim = 10)

res1
summary(res1)

## compare
f2 <- compare_sim_formulas("3-lvl" = f, "2-lvl" = f1)
res2 <- simulate(p, formula = f2, nsim = 50, cores = 1, satterthwaite = TRUE)
summary(res2)

## '2-lvl' should not give Power_bw = NA
# renome 'satterth_test' to test? also use for CI

p <- update(p, effect_size = cohend(1))


f1 <- sim_formula("y ~ treatment + pretest + (1 | cluster)", data_transform = transform_to_posttest, satterth_test = "treatment")
f2 <- compare_sim_formulas("3-lvl" = f, "2-lvl" = f1)
res2 <- simulate(p, formula = f2, nsim = 500, cores = 10, satterthwaite = TRUE)
summary(res2)

## TODO: fix so summary show all terms, e.g. 'pretest'

