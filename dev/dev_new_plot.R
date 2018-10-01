p <- study_parameters(n1 = 5, n2 = 1000, icc_pre_subject = 0.5, var_ratio = 0.03)

d <- simulate_data(p)

ggplot(d, aes(time, intercept_subject + slope_subject * time, group = subject, color = treatment)) + geom_line() +
    facet_grid(~treatment)


p <- study_parameters(n1 = 5, n2 = 1, n3 = 1e5, icc_pre_subject = 0.5, var_ratio = 0.03, icc_slope = 0.05)

d <- simulate_data(p)


ggplot(d, aes(time, intercept_cluster + slope_cluster * time, group = subject, color = treatment)) + geom_line() +
    facet_grid(~treatment)


d$mu <- with(d, intercept_cluster + slope_cluster * time)


tx <- lapply(unique(d$time), function(i, group) {
    x <- eta_sum(d[d$treatment == group & d$time == i, "mu"])
    x <- as.data.frame(x)
    x$treatment <- group
    x$time <- i

    x
}, group = 1)
tx <- do.call(rbind, tx)


cc <- lapply(unique(d$time), function(i, group) {
    x <- eta_sum(d[d$treatment == group & d$time == i, "mu"])
    x <- as.data.frame(x)
    x$treatment <- group
    x$time <- i

    x
}, group = 0)
cc <- do.call(rbind, cc)

x <- rbind(tx, cc)


tmp <- x[, !colnames(x) %in% c("var","mean", "sd", "Q50")]
tmp <- lapply(list(c("Q0.5", "Q99.5"),
                   c("Q2.5", "Q97.5"),
                   c("Q10", "Q90"),
                   c("Q25", "Q75")),
              function(x) {
                  d <- tmp[, c("treatment","time", x)]
                  data.frame(treatment = d$treatment,
                             time = d$time,
                             min = d[, x[1]],
                             max = d[, x[2]],
                             width = switch(x[1],
                                            "Q0.5" = 0.99,
                                            "Q2.5" = 0.95,
                                            "Q10" = 0.8,
                                            "Q25" = 0.5)
                  )
              })
tmp <- do.call(rbind, tmp)
tmp$width <- factor(tmp$width,
                    levels = rev(c(0.5, 0.8, 0.95, 0.99)))

ggplot(x, aes(time, mean, group = treatment)) +
    geom_ribbon(data = tmp, aes(ymin = min, ymax = max, y = NULL, x = time, group = width, fill = width), alpha = 0.75) +
    geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
    geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
    geom_point(aes(y = Q50), color = "red") +
    scale_color_manual(values = c("median" = "red", "mean" = "red")) +
    scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
    labs(linetype = "", color = "") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    facet_wrap(~treatment, ncol = 2) +
    scale_fill_brewer() +
    theme_minimal()
