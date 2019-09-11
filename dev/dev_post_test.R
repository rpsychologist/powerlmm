# Allow post process fitted object
# e.g. general linear hypothesis
# TODO: * function that creates return df

library(emmeans)
library(lmerTest)
p <- study_parameters(n1 = 3,
                      n2 = 20,
                      n3 = 4,
                      T_end = 10,
                      icc_pre_subject = 0.5,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

d <- simulate_data(p)
d$cluster <- factor(d$cluster)
d$treatment <- factor(d$treatment)
fit <- lmer(y ~ time*cluster + (1 + time | subject), data = d)


res <- emmeans(fit, ~cluster, at = list(time = 10), lmer.df = "satterthwaite")

contrast(res, list("cluster1" = c(rep(1/4, 4), rep(-1/4, 4)) ))

contrast(res)

d %>%
    filter(time == 10) %>%
    group_by(treatment) %>%
    summarise(y = mean(y)) %>%
    ungroup() %>%
    summarise(diff(y))


fit2 <- lmer(y ~ time*treatment + (1 + time | subject), data = d)




post_test <- function(fit, d = NULL) {
    res <- emmeans::emmeans(fit, ~cluster * time,
                   at = list(time = 10),
                   lmer.df = "satterthwaite")

    L <- c(rep(1/4, 4), rep(-1/4, 4))
    out <- emmeans::contrast(res,
                    list("TE" = L))
    out <- as.data.frame(out)

    out <-        data.frame(parameter = out$contrast,
                             estimate = out$estimate,
                             se = out$SE,
                             pval = out$p.value,
                             df = out$df,
                             df_bw = NA)


    # test slopes
    slope <- lmerTest::contest1D(fit, L = c(rep(0, 8), L))
    slope <- data.frame(parameter = "TE_slope_diff",
               estimate = slope$Estimate,
               se = slope$`Std. Error`,
               pval = slope$`Pr(>|t|)`,
               df = slope$df,
               df_bw = NA)

    rbind(out, slope)

}





f <- sim_formula("y ~ 0 + factor(cluster) + time:factor(cluster) + (1 + time | subject)",
                 test = NULL,
                 post_test = post_test)

p <- study_parameters(design = study_design(),
                      n1 = 3,
                      n2 = 20,
                      n3 = 4,
                      T_end = 10,
                      icc_pre_subject = 0.5,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

res <- simulate(p, nsim = 1000, formula = f, cores = 16)

summary(res)
get_power(p)

# Fixed cluster DGP
fixed_cluster <- function(d) {
    d$cluster_slope <- c(0.5, 0.25, -0.25, -0.5,
                         0.75, 0.5, -0.5, -0.75)[d$cluster]
    # d$cluster_slope <- c(0.5, 0.25, -0.25, -0.5,
    #                      0.75, 0.5, 0.2, 0.1)[d$cluster]

    d$y <- d$y + d$time * d$cluster_slope
    d
}


simulate_data(p) %>%
    fixed_cluster %>%
    group_by(time, cluster, treatment) %>%
    summarise(y = mean(y)) %>%
    ggplot(aes(time, y, group = cluster, color = cluster, linetype = factor(treatment))) + geom_line()




f <- sim_formula("y ~ 0 + factor(cluster) + time:factor(cluster) + (1 + time | subject)",
                 test = NULL,
                 data_transform = fixed_cluster,
                 post_test = post_test)

f2 <- sim_formula("y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)",
                  test = NULL,
                  data_transform = fixed_cluster)

f3 <- sim_formula("y ~ time*treatment + (1 + time | subject)",
                  test = NULL,
                  data_transform = fixed_cluster)

p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_slope = 0,
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = 0)

res <- simulate(p, nsim = 1000,
                formula = sim_formula_compare("fixed" = f,
                                              "random" = f2,
                                              "ignore" =f3), cores = 16)

summary(res)
