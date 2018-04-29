library(lme4)

p <- study_parameters(n1 = 11,
                      n2 = 20,
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


f0 <- sim_formula("y ~ treatment",
            test = "treatment",
            data_transform = transform_to_posttest)

f1 <- sim_formula("y ~ treatment + (1 | cluster)",
                  test = "treatment",
                  data_transform = transform_to_posttest)

res <- simulate(p, nsim = 50, formula = sim_formula_compare("m0" = f0, "m1" = f1),
                cores = 1,
                satterthwaite = TRUE)

summary(res)
summary(res, model_selection = "FW")

summary(res, df_bw = list("treatment" = 4))
summary(res, model_selection = "FW")



res1 <- simulate(p, nsim = 50, formula = sim_formula("y ~ treatment + (1 | cluster)",
                                                   test = "treatment",
                                                   data_transform = transform_to_posttest),
                 CI = TRUE)
summary(res1)




res1 <- simulate(p, nsim = 50, formula = sim_formula("y ~ time*treatment + (1 | cluster)",
                                                     test = c("treatment", "time", "time:treatment")),
                 satterthwaite = TRUE,
                 CI = TRUE)
summary(res1)
summary(res1, df_bw = list("treatment" = 20, "time" = 20, "time:treatment" = 20))


## Multi

p <- update(p, n2 = 20:21)

f0 <- sim_formula("y ~ treatment",
                  test = "treatment",
                  data_transform = transform_to_posttest)

f1 <- sim_formula("y ~ treatment + (1 | cluster)",
                  test = "treatment",
                  data_transform = transform_to_posttest)

res <- simulate(p, nsim = 50, formula = sim_formula_compare("m0" = f0, "m1" = f1),
                cores = 1,
                satterthwaite = TRUE,
                CI = FALSE)

summary(res, para = "treatment", model = "m0")
summary(res, para = "treatment", model_selection = "FW")
summary(res, para = "treatment", model = "m0", df_bw = list("treatment" = 1))
x <- summary(res, model_selection = "FW", para = "treatment")
