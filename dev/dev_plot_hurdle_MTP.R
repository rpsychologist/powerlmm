p_mtp <- study_parameters(
    design = structure(list(), class = "plcp_hurdle"),
    n1 = 3,
    n2 = 20,
    fixed_intercept = log(300), # median(Y > 0)
    fixed_hu_intercept = qlogis(0.5), # prop == 0
    fixed_slope = log(0.99),
    fixed_hu_slope = log(1),
    sd_hu_intercept = 0.5,
    sd_hu_slope = 0,
    sd_intercept = 3,
    sd_slope = 0.1,
    cor_intercept_slope = -0.15,
    cor_intercept_hu_intercept = -0.66,
    cor_intercept_hu_slope = 0.2,
    cor_slope_hu_intercept = -0.1,
    cor_slope_hu_slope = -0.1,
    cor_hu_intercept_hu_slope = 0.15,
    shape = 1.6,
    RR_cont = 0.8,
    OR_hu = 2,
    marginal = TRUE,
    family = "gamma")

m <- marginalize(p_mtp, R = 1e5)


# time
## overall
plot_hurdle_time(m$y_overall) + scale_y_log10()
plot_hurdle_time(m$y_positive) + scale_y_log10()
plot_hurdle_time(m$hu_prob)


# Diff
plot_hurdle_diff(m)
plot_hurdle_diff(m, fixed_overall = get_overall_hurdle(p_mtp))

plot_hurdle_diff(m, hu = TRUE)

plot_hurdle_probs(m1)



tmp <- m1$mu_overall_vec
ps <- seq(0.00000000001, 0.99999999, length.out = 100)
y <- unlist(lapply(ps, function(p) quantile(tmp[[6]], p)/quantile(tmp[[3]], p)))


plot(ps, y)

median(tmp[[6]])/median(tmp[[3]])

mean(tmp[[6]]) - mean(tmp[[3]])
mean(tmp[[6]]) / mean(tmp[[3]])



quantile(tmp[[6]], 0.78)/quantile(tmp[[3]], 0.78)

plot(density(tmp[[6]]/tmp[[3]]))
