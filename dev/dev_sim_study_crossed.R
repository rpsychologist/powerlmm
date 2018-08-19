res <- simulate(p,
                formula = sim_formula("y ~ time*treatment + (time | subject) + (time * treatment | cluster)"),
                nsim = 5)

res
summary(res)



d <- simulate_data(p)
fit <- lmer(y ~ time*treatment + (time | subject) + (time * treatment | cluster), data = d)

as.data.frame(VarCorr(fit))
