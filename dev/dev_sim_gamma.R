library(tidyverse)
library(brms)
des <- study_design(family = "gamma")

p <- study_parameters(design = des,
                      n1 = 11,
                      n2 = 50,
                      fixed_intercept = 5,
                      sigma_subject_intercept = 1,
                      effect_size = log(1),
                      shape = 2)

d <- simulate_data(p)

d$p <- exp(d$fixed_intercept + d$subject_intercept + (d$subject_slope + d$fixed_slope) * d$time)

d %>% group_by(time, treatment) %>%
    summarise(mean(y))


ggplot(d, aes(time, p, group = subject, color = treatment)) + geom_line()

fit_brm <- brm(y ~ time * treatment + (1 | subject),
               data = d,
               iter = 1,
               chains = 1,
               family = Gamma("log"))

f_brm <- sim_formula(fit_brm, iter = 1000, silent = TRUE)

f_glmer <- sim_formula("y ~ time * treatment + (1 | subject)", family = Gamma("log"))
# f_glmer <- sim_formula(create_lmer_formula(p))


res <- simulate(p,
                formula = sim_formula_compare("brms" = f_brm,
                                              "glmer" = f_glmer),
                nsim = 100,
                cores = 16)


summary(res)
summary(res, model = "brms")


#
f_glmer <- sim_formula("y ~ time * treatment + (1 | subject)", family = Gamma("log"), nAGQ=0)

res2 <- simulate(p,
                formula = sim_formula_compare("glmer" = f_glmer),
                nsim = 100,
                cores = 16)


summary(res2)
