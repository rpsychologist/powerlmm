library(lme4)

p <- study_parameters(n1 = 11,
                      n2 = 4,
                      n3 = 4,
                      icc_pre_cluster = 0.1,
                      icc_pre_subject = 0.5,
                      )

d <- simulate_data(p)
d <- transform_to_posttest(d)

fit0 <- lm(y ~ treatment, d)
fit1 <- lme4::lmer(y ~ treatment + (1 | cluster), data = d)

anova(fit1, fit0, refit = FALSE)
ranova(fit1)

logLik(fit0, REML = TRUE)
logLik(fit1, REML = TRUE)

##
res <- simulate(p, nsim = 50, formula = sim_formula("y ~ treatment",
                                                   test = "treatment",
                                         data_transform = transform_to_posttest),
                satterthwaite = TRUE)

summary(res)

res1 <- simulate(p, nsim = 50, formula = sim_formula("y ~ treatment + (1 | cluster)",
                                                   test = "treatment",
                                                   data_transform = transform_to_posttest),
                 satterthwaite = TRUE)
summary(res1)

