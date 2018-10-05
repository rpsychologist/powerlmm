# Shows how to use data_transform to generate MAR missingness
# Compares LMM vs GEE

p <- study_parameters(n1 = 11,
                      n2 = 100,
                      icc_pre_subject = 0.6,
                      fixed_slope = -2,
                      var_ratio = 0.02,
                      cor_subject = 0,
                      effect_size = 0)




add_MAR_missing <- function(data) {
    d <- data

    time <- unique(d$time)
    for(i in seq_along(time)) {
        for(j in unique(d$treatment)) {
            if(i == 1) {
                p_miss <- 0
            } else {
                y <- d$y[d$time == time[i-1] & d$treatment == j]
                y <- y - mean(y)
                p_miss <- pnorm(-2 + -0.05*y + 0.01*time[i] + 0.1 *j *y)
            }

            d[d$time == time[i] & d$treatment == j, "p_miss"] <- p_miss
        }
    }

    d$miss <- rbinom(nrow(d), 1, d$p_miss)

    for(i in unique(d$subject)) {
        tmp <- d[d$subject == i, ]
        dropout <- which(tmp$miss == 1)[1]
        if(!is.na(dropout)) tmp[dropout:nrow(tmp), "y"] <- NA

        d[d$subject == i, "y"] <- tmp$y
    }
    d
}
d <- simulate_data(p)
d <- add_MAR_missing(d)

d %>%
    group_by(treatment, time) %>%
    summarise(mean(is.na(y))) %>%
    print(n=50)

d %>% filter(time == 10) %>%
    group_by(treatment) %>%
    summarise(mean(y, na.rm=TRUE),
              mean(y_c))


 ggplot(d, aes(time, p_miss, group = subject)) + geom_line() +
    facet_grid(~treatment)


 d$pre <- rep(d$y[d$time == 0], each = 11)
 d %>%
     arrange(pre) %>%
     mutate(subject = factor(subject, levels = unique(subject))) %>%
 ggplot(aes(time, subject, fill = is.na(y))) + geom_tile() +
     facet_wrap(~treatment, scales = "free", drop = TRUE)

 p2 <- p
 p2$n2 <- 2000
 d <- simulate_data(p2)
 d <- add_MAR_missing(d)
 d$mu <- p$fixed_intercept + d$subject_intercept + (p$fixed_slope + d$subject_slope) * d$time
 #d$mu[is.na(d$y)] <- NA
 ggplot(d, aes(time, mu, group = subject, color = is.na(y))) +
     geom_line(aes(alpha = is.na(y))) +
     stat_summary(geom = "line", fun.data.y= mean, aes(y = y, group = treatment), linetype = "dotted", color = "blue", size = 1) +
     stat_summary(geom = "line", fun.data.y= mean, aes(y = y_c, group = treatment), color = "blue", size = 1) +
     facet_grid(~treatment) +
     scale_alpha_manual(values = c(0.33, 0)) +
     scale_color_manual(values = c("black", "red")) +
     theme_minimal()


# Simulate ----------------------------------------------------------------
f <- sim_formula("y ~ time*treatment + (1 | subject)", data_transform = add_MAR_missing, family = gaussian())
f0 <- sim_formula("y ~ time*treatment + (1 | subject)", family = gaussian())

 res <- simulate(p, formula = sim_formula_compare("MAR" = f,
                                                  "complete" = f0), nsim = 1000, cores = 4)
summary(res)

# GEE
library(geepack)
library(future.apply)
plan(multiprocess, workers = 4)
sim_gee <- function(i) {
    d <- simulate_data(p)
    d <- add_MAR_missing(d)
    fit <- geeglm(y ~ time*treatment, id = subject, data = d)
    x <- summary(fit)
    x$coefficients[4,]
}
res_gee <- future_lapply(1:1000, sim_gee)
res_gee <- do.call(rbind, res_gee)

mean(res_gee$Estimate)
mean(res_gee$`Pr(>|W|)` < 0.05)


# Obs mean
sim_obs_mean_diff <- function(i) {
    d <- simulate_data(p)
    d <- add_MAR_missing(d)
    d <- d[d$time == max(d$time), ]
    tx <- d[d$treatment == 1, ]
    cc <- d[d$treatment == 0, ]

    data.frame("MAR" = mean(tx$y, na.rm=TRUE) - mean(cc$y, na.rm=TRUE),
         "complete" = mean(tx$y_c) - mean(cc$y_c),
         "miss_tx" = mean(is.na(tx$y)),
         "miss_cc" = mean(is.na(cc$y)))

}
res_obs_diff <- future_lapply(1:1000, sim_obs_mean_diff)
res_obs_diff <- do.call(rbind, res_obs_diff)
colMeans(res_obs_diff)
