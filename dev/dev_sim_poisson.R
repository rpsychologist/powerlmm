library(tidyverse)
des <- study_design(family = "poisson")

p <- study_parameters(design = des,
                      n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      fixed_intercept = log(5),
                      fixed_slope = log(0.95),
                      var_ratio = 0.02,
                      effect_size = log(1),
                      sigma_error = 1)

d <- simulate_data(p)
d$p <- exp(d$fixed_intercept + d$subject_intercept + (d$subject_slope + d$fixed_slope) * d$time)

d %>% group_by(time, treatment) %>%
    summarise(mean(y))


ggplot(d, aes(time, p, group = subject, color = treatment)) + geom_line()

res <- simulate(p, nsim = 500, cores = 16)

summary(res)

res2 <- simulate(p, nsim = 16, cores = 16, CI = TRUE)



# fit LMM instead
f <- sim_formula("y ~ time*treatment + (1 + time | subject)")
res3 <- simulate(p, formula = f, nsim = 500, cores = 16)
summary(res3)


# gaussian
p <- study_parameters(n1 = 11,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      effect_size = log(2),
                      sigma_error = 1)
p$family
