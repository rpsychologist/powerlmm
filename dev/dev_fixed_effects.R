









# New argument for fixed clusters DGP -------------------------------------
## want to use 3-level RE model, but be able to set fixed cluster based on percentiles AND sds?
## create helper that specify fixed cluster slopes in cohen's d units?
##          - just use the existing cohens func.


# treatment effect will be unaffected as long as the therapist skill is same in both groups
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      fixed_slope = -1,
                      icc_slope = 0.1,
                      fixed_cluster_intercepts = list(c(0.5,   0.5,   0.5,   0.5)),
                      fixed_cluster_slopes     = list(c(0.9, 0.9, 0.9, 0.9)), # equal in both groups
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = -1)


d <- simulate_data(p)
d$cluster_slope %>% unique

# use per treatment to set different therapist slopes in tx and cc

p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_slope = 0.1,
                      fixed_cluster_intercepts = 0.5, # default to 0
                      fixed_cluster_slopes = per_treatment(c(0.8, 0.85, 0.9, 0.95),
                                                           c(0.1, 0.2, 0.25, 0.3)),
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

d <- simulate_data(p)
d$cluster_slope %>% unique
d$cluster_intercept %>% unique

## per_treatment single
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0.1,
                      icc_slope = 0.1,
                      fixed_cluster_intercepts = per_treatment(0.5, 0.1), # default to 0
                      fixed_cluster_slopes = per_treatment(c(0.8, 0.85, 0.9, 0.95),
                                                           c(0.1, 0.2, 0.25, 0.3)),
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

d <- simulate_data(p)
d$cluster_intercept %>% unique
d$cluster_slope %>% unique


## fixed effects

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

