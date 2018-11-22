## fixed effects
## TODO: pass L via sim_formula(test = L)?


library(emmeans)
library(lmerTest)
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

d <- simulate_data(p)
d <- trans_f2(d)
d$cluster <- factor(d$cluster)

fit0 <- lm(y ~ treatment, data = d)
summary(fit0)

fit <- lm(y ~ 0 + cluster, data = d)
summary(fit)



fm <- lmer(y ~ treatment + (1 + treatment | cluster), data = d)
summary(fm)

K <- diag(length(coef(fit)))[-1,]
rownames(K) <- names(coef(fit))[-1]
K

summary(glht(fit, linfct = rbind("delta" = c(1/4, 1/4, 1/4, 1/4, -1))))


d %>% group_by(cluster, treatment) %>%
    summarise(y = mean(y)) %>%
    group_by(treatment) %>%
    summarise(mean(y))

d %>% group_by(treatment) %>%
    summarise(mean(y))

# lmer
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      partially_nested = TRUE,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      cor_subject = -0.5,
                      effect_size = 0)

d <- simulate_data(p)

d[d$treatment == 0, ]$cluster <- 1
d$cluster <- as.factor(d$cluster)

# with random slope
fit0 <- lmer(y ~ 0 + cluster + cluster:time + (1 + time | subject), data = d)
summary(fit0)

contest1D(fit0, L = c(0,0,0,0,0, -1, 1/4, 1/4, 1/4, 1/4))

# without random slope
fit1 <- lmer(y ~ 0 + cluster + cluster:time + (1 | subject), data = d)

contest1D(fit1, L = c(0,0,0,0,0, -1, 1/4, 1/4, 1/4, 1/4))


fit2 <- lmer(y ~ 0 + treatment * time + (1 + time | subject), data = d)
summary(fit2)

