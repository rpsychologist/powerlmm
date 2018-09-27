library(ggplot2)

des <- structure(list(), class = "plcp_hurdle")

p1 <- study_parameters(
    design = des,
    n1 = 3,
    n2 = 20,
    fixed_intercept = log(30), # median(Y > 0)
    fixed_hu_intercept = qlogis(0.8), # prop == 0
    fixed_slope = log(0.99),
    fixed_hu_slope = log(1),
    sd_hu_intercept = 3,
    sd_hu_slope = 0.2,
    sd_intercept = 2,
    sd_slope = 0.05,
    cor_intercept_slope = -0.15,
    cor_intercept_hu_intercept = -0.66,
    cor_intercept_hu_slope = 0.2,
    cor_slope_hu_intercept = -0.1,
    cor_slope_hu_slope = -0.1,
    cor_hu_intercept_hu_slope = 0.15,
    shape = 1.6,
    RR_cont = 0.8,
    OR_hu = 2,
    marginal = FALSE,
    family = "gamma")


m1 <- marginalize(p1, R = 1e5)


# time
## overall
plot_hurdle_time(m1$y_overall)
plot_hurdle_time(m1$y_positive) + scale_y_log10()
plot_hurdle_time(m1$hu_prob)


# Diff
plot_hurdle_diff(m1)
plot_hurdle_diff(m1, hu = TRUE)

plot_hurdle_probs(m1)



ggplot(tmp, aes(percentile, p, color = percentile == 0.5)) +
    geom_point() +
    geom_segment(aes(y = 0, yend = OR, xend = percentile)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = m$post_hu[m$post_hu$var == "marg_OR", "est"]) +
    theme_minimal()



#
x <- seq(-3, 3, length.out = 100)
plot(x, plogis(x+1) - plogis(x+1))
