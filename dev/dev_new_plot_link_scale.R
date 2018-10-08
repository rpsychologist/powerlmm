p_mtp <- study_parameters(
    design = study_design(family = "hurdle"),
    n1 = 11,
    n2 = 1e5,
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


d <- simulate_data(p_mtp)

d$mu2 <- with(d, subject_intercept + subject_slope * time)

# 3-level
d$mu3 <- with(d, cluster_intercept + cluster_slope * time)

# lvl 2
cc <- eta_sum_d(d = d, var = "mu2", treatment = 0)
tx <- eta_sum_d(d = d, var = "mu2", treatment = 1)
x <- rbind(cc, tx)
Q_long <- reshape_eta_sum(x)

# lvl 3
cc3 <- eta_sum_d(d = d, var = "mu3", treatment = 0)
tx3 <- eta_sum_d(d = d, var = "mu3", treatment = 1)
x3 <- rbind(cc3, tx3)
Q_long3 <- reshape_eta_sum(x3)

ymin <- min(Q_long$min, Q_long3$min)
ymax <- max(Q_long$max, Q_long3$max)

p2 <- ggplot(x, aes(time, mean, group = treatment)) +
    geom_ribbon(data = Q_long, aes(ymin = min, ymax = max, y = NULL, x = time, group = interaction(width, treatment), fill = width), alpha = 0.75) +
    geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
    geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
    geom_point(aes(y = Q50), color = "red") +
    scale_color_manual(values = c("median" = "red", "mean" = "red")) +
    scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
    labs(linetype = "", color = "") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    facet_wrap(~treatment, ncol = 2) +
    lims(y = c(ymin, ymax)) +
    scale_fill_brewer() +
    theme_minimal()




p3 <- ggplot(x, aes(time, mean, group = treatment)) +
    geom_ribbon(data = Q_long3, aes(ymin = min, ymax = max, y = NULL, x = time, group = interaction(width, treatment), fill = width), alpha = 0.75) +
    geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
    geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
    geom_point(aes(y = Q50), color = "red") +
    scale_color_manual(values = c("median" = "red", "mean" = "red")) +
    scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
    labs(linetype = "", color = "") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    facet_wrap(~treatment, ncol = 2) +
    lims(y = c(ymin, ymax)) +
    scale_fill_brewer() +
    theme_minimal()

gridExtra::grid.arrange(p2,p3)
