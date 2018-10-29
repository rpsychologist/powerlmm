

geom_path(data = trend$x,
          aes(x = mean,
              y = time,
              linetype = "mean",
              fill = NULL,
              group = interaction(treatment, var)),
          size = 1) +

plot(m, RE = TRUE, type = "trend_ridges",
     RE_level = c(1,2), trim = c(0, 1),
     sd2_p = c(0.5), sd2_hu_p = c(0.5)) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 100, 1000, 10000)) +

    scale_color_manual(values = c("#222f3e", "#ff6b6b"))



p <- study_parameters(
    design = study_design(family = "hurdle"),
    n1 = 3,
    n2 = 20,
    fixed_intercept = log(100), # median(Y > 0)
    fixed_hu_intercept = qlogis(0.2), # prop == 0
    fixed_slope = log(0.99),
    fixed_hu_slope = log(1),
    sd_hu_intercept = 1,
    sd_hu_slope = 0.2,
    sd_intercept = 0.3,
    sd_slope = 0.05,
    cor_intercept_slope = -0.15,
    cor_intercept_hu_intercept = -0.66,
    cor_intercept_hu_slope = 0.2,
    cor_slope_hu_intercept = -0.1,
    cor_slope_hu_slope = -0.1,
    cor_hu_intercept_hu_slope = 0.15,
    shape = 2,
    RR_cont = 0.33,
    OR_hu = 1,
    marginal = TRUE,
    family = "gamma")

m <- marginalize(p)


res <- .mu_vec_to_long(m,
                       RE_level = 1,
                       var1 = "mu1_vec",
                       var2 = "mu_overall_vec",
                       tx_var = "y_overall",
                       level1_func = .sample_level1_nested_hurdle, R = 1e5)

tmp <- dplyr::filter(res, time == 0)
y <- compute_density(filter(tmp, y > 0)$y)
y$density <- y$count/nrow(tmp)
y <- subset(y, x > 0)

y0 <- data.frame(x = 0,
                 density = 0,
                 count = length(dplyr::filter(res, time == 0, y - 0 < .Machine$double.eps^0.5)$y))
y0$density <- y0$count/nrow(tmp)

#y  <- rbind(y, y0)

ggplot(y, aes(x, density)) +
    geom_density(stat = "identity", fill = "black") +
    geom_segment(data = y0, aes(y = 0, yend = density, x = 0, xend = 0), stat = "identity", fill = "black", size = 2)
   # scale_x_continuous(trans = "log1p") +
   # scale_y_continuous(trans = "sqrt")

