p <- study_parameters(n1 = 5, n2 = 1e5, icc_pre_subject = 0.5, var_ratio = 0.03, cor_subject = -0.8)

d <- simulate_data(p)

ggplot(d, aes(time, intercept_subject + slope_subject * time, group = subject, color = treatment)) + geom_line() +
    facet_grid(~treatment)


p <- study_parameters(n1 = 5, n2 = 1, n3 = 1e5,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.03,
                      icc_slope = 0.05,
                      icc_pre_cluster = 0.1,
                      partially_nested = TRUE,
                      effect_size = cohend(0.5))

d <- simulate_data(p)


ggplot(d, aes(time, intercept_cluster + slope_cluster * time, group = subject, color = treatment)) + geom_line() +
    facet_grid(~treatment)


# 2-level
d$mu <- with(d, intercept_subject + slope_subject * time)


# 3-level
d$mu3 <- with(d, intercept_cluster + slope_cluster * time)


cc <- eta_sum_d(d = d, var = "mu", treatment = 0)
tx <- eta_sum_d(d = d, var = "mu", treatment = 1)
x <- rbind(cc, tx)
Q_long <- reshape_eta_sum(x)

ggplot(x, aes(time, mean, group = treatment)) +
    geom_ribbon(data = Q_long, aes(ymin = min, ymax = max, y = NULL, x = time, group = interaction(width, treatment), fill = width), alpha = 0.75) +
    geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
    geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
    geom_point(aes(y = Q50), color = "red") +
    scale_color_manual(values = c("median" = "red", "mean" = "red")) +
    scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
    labs(linetype = "", color = "") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    #facet_wrap(~treatment, ncol = 2) +
    scale_fill_brewer() +
    theme_minimal()



#
ggplot(x, aes(time, mean, group = treatment)) +
    geom_ribbon(data = Q_long, aes(ymin = min, ymax = max, y = NULL, x = time, group = interaction(width, treatment), fill = width), alpha = 0.5) +

    geom_line(data = Q_long, aes(y = max, x = time, group = interaction(width, treatment), alpha = width), fill = NA, color = "#08519C") +
    geom_line(data = Q_long, aes(y = min, x = time, group = interaction(width, treatment), alpha = width),  fill = NA, color = "#08519C") +
    geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
    geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
    geom_point(aes(y = Q50), color = "red") +
    scale_color_manual(values = c("median" = "red", "mean" = "red")) +
    scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
    labs(linetype = "", color = "") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    #facet_wrap(~treatment, ncol = 2) +
    scale_fill_brewer() +
    theme_minimal()


# percentiles

ps <- 1:99/100

mu <- d[d$treatment == 0 & d$time == max(d$time), "mu"]
post_cc <- data.frame("percentile" = ps,
                   "value" = quantile(mu, ps),
                   "treatment" = 0
)

mu <- d[d$treatment == 1 & d$time == max(d$time), "mu"]
post_tx <- data.frame("percentile" = ps,
                      "value" = quantile(mu, ps),
                      "treatment" = 1
)

post_tx$diff <- post_tx$value - post_cc$value

ES <- get_slope_diff(p)

ggplot(post_tx, aes(percentile, diff)) +
    geom_histogram(stat = "identity", color = "white", fill = "#3498db", alpha = .75) +
    ylim(c(-10, 10)) +
    geom_hline(yintercept = ES, linetype = "dotted", alpha = 0.75, size = 0.75)



    geom_hline(yintercept = ES_med, linetype = "dashed", alpha = 0.75, size = 0.75) +
    scale_y_continuous(sec.axis = sec_axis(~ ., breaks = c(ES_med, ES),
                                           labels = c(paste(round(ES_med, 2), " (median)"),
                                                      paste(round(ES, 2), " (mean)")
                                           )
    )
    ) +        theme_minimal() +
    theme(legend.position = "none")
