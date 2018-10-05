des <- study_design(family = "binomial")

p <- study_parameters(design = des,
                      n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.02,
                      effect_size = log(1),
                      sigma_error = 1)

d <- simulate_data(p)
d$p <- plogis(d$fixed_intercept + d$subject_intercept + (d$subject_slope + d$fixed_slope) * d$time)

ggplot(d, aes(time, p, group = subject, color = treatment)) + geom_line()

res <- simulate(p, nsim = 500, cores = 16)

summary(res)



# multi -------------------------------------------------------------------

des <- study_design(family = "binomial")

p <- study_parameters(design = des,
                      n1 = c(5,11),
                      n2 = 25,
                      T_end = 10,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.02,
                      effect_size = log(1),
                      sigma_error = 1)

f <- sim_formula("y ~ time*treatment + (1 + time | subject)", family = binomial("logit"))

res <- simulate(p, formula = f, nsim = 16, cores = 16)
x <- summary(res)
