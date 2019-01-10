
# TODO

# * use correct labels based on family type

# Gaussian ----------------------------------------------------------------
p <- study_parameters(design = study_design(),
                          n1 = 11,
                          n2 = 25,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.02,
                      cor_subject = -0.5,
                          effect_size = log(0.5),
                          sigma_error = 1)

m <- marginalize(p)

plot(m)
plot(p , RE = TRUE, type = "trend", RE_level = c(1,2))
plot(p , RE = TRUE, type = "trend", RE_level = c(2,3))

plot(m , RE = TRUE, type = "trend", RE_level = c(1,2))


# Set better limits?
# trim data to 99%?
plot(m, type = "trend_ridges", RE_level = c(1, 2,3)) + facet_wrap(~var, scales = "free_x")

# TODO: support link scale
plot(p, type = "trend_ridges", RE_level = c(1, 2,3))

plot(m , RE = FALSE, type = "post_diff", RE_level = c(2))
plot(m , RE = FALSE, type = "post_diff", RE_level = c(3))

# Binomial ----------------------------------------------------------------
p_bin <- study_parameters(design = study_design(family = "binomial"),
                          n1 = 11,
                          n2 = 25,
                          fixed_intercept = qlogis(0.7),
                          sigma_subject_intercept = sqrt(pi^2/3),
                          sigma_subject_slope = 0,
                          cor_subject = -0.5,
                          effect_size = log(0.5),
                          sigma_error = sqrt(pi^2/3))

## link scale
plot(p_bin, RE_level = c(1, 2))

plot(p_bin,
     RE = FALSE,
     type = "trend",
     RE_level = c(2, 3))


plot(p_bin,
     RE = TRUE,
     type = "trend",
     RE_level = c(2, 3))

plot(p_bin,
     RE = FALSE,
     type = "trend_ridges",
     RE_level = c(2, 3))


# show warning
plot(p_bin, type = "post_ratio")
plot(p_bin, type = "post_diff")


## marginalized
m_bin <- marginalize(p_bin)

plot(m_bin, RE_level = c(1, 2))

## TODO: fix level 1
plot(m_bin,
     RE = FALSE,
     type = "trend",
     RE_level = c(1,2,3))

plot(m_bin,
     RE = TRUE,
     type = "trend",
     RE_level = c(1,2,3))

plot(m_bin,
     RE = FALSE,
     type = "trend_ridges",
     RE_level = c(1,2,3)) +
    facet_wrap(~var, nrow  = 3)

plot(m_bin, type = "post_ratio")
plot(m_bin, type = "post_diff")


plot(p_bin, type = "post_ratio")
plot(p_bin, type = "post_diff")

# Poisson ----------------------------------------------------------------
p_pois <- study_parameters(design = study_design(family = "poisson"),
                          n1 = 4,
                          n2 = 25,
                          T_end = 10,
                          fixed_intercept = log(10),
                          sigma_subject_intercept = 1,
                          sigma_subject_slope = 0.02,
                          cor_subject = -0.5,
                          sigma_cluster_intercept = 0.4,
                          effect_size = log(0.5),
                          sigma_error = 1)

m_pois <- marginalize(p_pois)
plot(p_pois)
plot(m_pois)
plot(m_pois, type = "post_ratio")

plot(m_pois, type = "trend", RE_level = c(1,2))

# TODO:
# * fix lims
plot(m_pois, type = "trend_ridges", RE_level = c(1,2)) + xlim(0, 20)

# lognormal ----------------------------------------------------------------
p_ln <- study_parameters(design = study_design(family = "lognormal"),
                           n1 = 11,
                           n2 = 25,
                           fixed_intercept = log(1000),
                           sigma_subject_intercept = 1,
                           sigma_subject_slope = 0.02,
                         sigma_cluster_intercept = 0.4,
                           effect_size = log(1),
                           sigma_error = 1)

m_ln <- marginalize(p_ln)
plot(m_ln)
plot_link(p_ln)
plot(m_ln, type = "trend", RE_level = c(1,2))
plot(m_ln, type = "trend_ridges", RE_level = c(1,2)) + xlim(0, 10000) + scale_x_log10()

# Gamma ----------------------------------------------------------------
p_gamma <- study_parameters(design = study_design(family = "gamma"),
                         n1 = 3,
                         n2 = 25,
                         T_end = 10,
                         fixed_intercept = log(500),
                         sigma_subject_intercept = 1,
                         sigma_cluster_intercept = 0.2,
                         effect_size = log(0.5),
                         shape = 2,
                         sigma_error = 1)

m_gamma <- marginalize(p_gamma)

plot(p_gamma, RE_level = c(1,2,3))
plot(m_gamma)

plot(m_gamma, type = "trend", RE_level = c(1,2), sd2_p = c(0.5, 0.5)) +
    ylim(0, 5000) +
    scale_y_log10() +
        scale_color_brewer(palette = "Accent")

plot(m_gamma, type = "trend_ridges", RE_level = c(1,2)) +
    xlim(0, 3000) +
    scale_x_log10()


plot(m_gamma, type = "trend_ridges", RE_level = c(1,2,3))


# Hurdle models



sd(Intercept)                   1.84      0.29     1.31     2.46        709 1.00
sd(week)                        0.06      0.04     0.00     0.17        499 1.00
sd(hu_Intercept)                3.60      0.60     2.56     4.90       1944 1.00
sd(hu_week)                     0.29      0.09     0.15     0.49       1282 1.00
cor(Intercept,week)            -0.08      0.38    -0.81     0.66       2342 1.00
cor(Intercept,hu_Intercept)    -0.66      0.13    -0.86    -0.36        487 1.01
cor(week,hu_Intercept)         -0.16      0.41    -0.84     0.68        141 1.04
cor(Intercept,hu_week)          0.27      0.24    -0.21     0.71       2180 1.00
cor(week,hu_week)              -0.10      0.40    -0.79     0.71        365 1.02
cor(hu_Intercept,hu_week)       0.13      0.28    -0.46     0.63       2565 1.00

Population-Level Effects:
    Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
Intercept                           2.81      0.51     1.76     3.78        660 1.01
hu_Intercept                        2.72      0.68     1.48     4.15       1755 1.00
pre                                 0.83      0.43    -0.03     1.67       1982 1.00
treatmentBCTMbehandling             0.42      0.53    -0.59     1.44        962 1.00
week                               -0.01      0.05    -0.12     0.09        290 1.02
treatmentBCTMbehandling:week       -0.07      0.05    -0.19     0.02       1938 1.00
hu_pre                              0.92      0.62    -0.19     2.23       3020 1.00
hu_treatmentBCTMbehandling         -0.87      0.84    -2.53     0.75       1686 1.00
hu_week                             0.07      0.09    -0.10     0.28       1387 1.00
hu_treatmentBCTMbehandling:week     0.19      0.11    -0.01     0.41       2824 1.00
# CTP hurdle lognormal ----------------------------------------------------
library(cowplot)
p <- study_parameters(
    design = study_design(family = "hurdle"),
    n1 = 3,
    n2 = 20,
    T_end = 10,
    fixed_intercept = 2.81, # median(Y > 0)
    fixed_hu_intercept = qlogis(0.8), # prop == 0
    fixed_slope = log(1.3)/10,
    fixed_hu_slope = log(.5)/10,
    sd_hu_intercept = 3.6,
    sd_hu_slope = 0.3,
    sd_intercept = 1.84,
    sd_slope = 0.05,
    cor_intercept_slope = -0.1,
    cor_intercept_hu_intercept = -0.8,
    cor_intercept_hu_slope = 0.2,
    cor_slope_hu_intercept = -0.15,
    cor_slope_hu_slope = -0.1,
    cor_hu_intercept_hu_slope = 0.15,
    shape = 1.6,
    RR_cont = 0.5,
    OR_hu = 1.5,
    marginal = TRUE,
    family = "gamma")


# TODO: plot on linear scale
#plot(p)

m <- marginalize(p, R = 1e4)
RE <- FALSE
p0 <- plot(m, RE = RE, type = "trend",
     outcome = "hurdle",
     RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") + facet_wrap(~var)

p1 <- plot(m, RE = RE, type = "trend",
           outcome = "overall",
           RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 50, 500, 1000)) + facet_wrap(~var)


p2 <- plot(m, RE = RE, type = "trend",
           outcome = "positive",
           RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 50, 500, 1000, 1e4, 5e4)) + facet_wrap(~var)


plot_grid(p0,
          p1,
          p2,
          ncol = 1)

# RE EFFECTs
RE <- TRUE
p0 <- plot(m, RE = RE, type = "trend",
           outcome = "hurdle",
           RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") + facet_wrap(treatment~var)

p1 <- plot(m, RE = RE, type = "trend",
           outcome = "overall",
           RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 50, 500, 1000)) + facet_wrap(treatment~var)


p2 <- plot(m, RE = RE, type = "trend",
           outcome = "positive",
           RE_level = c(1,2)) +
    scale_fill_brewer(palette = "PuBu") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 50, 500, 1000, 1e4, 5e4)) + facet_wrap(treatment~var)


plot_grid(p0,
          p1,
          p2,
          ncol = 1)




#TODO remove tidyverse dependency
plot(m,
     RE = TRUE,
     type = "trend_ridges",
     RE_level = c(1,2), trim = c(0, 0.99),
     sd2_p = c(0.5, 0.5),
     sd2_hu_p = c(0.5)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 100, 500, 1000, 10000))


# Done
plot(m, RE = FALSE)
plot(m, RE = TRUE)

plot(m, RE = TRUE, type = "trend", RE_level = c(1,2)) + scale_fill_brewer(palette = "PuBu")
plot(m, RE = FALSE, type = "trend")

plot(m, RE = TRUE, type = "trend_ridges", RE_level = c(1, 2), trim = c(0, 1))

plot(m, RE = TRUE, type = "trend_ridges", RE_level = c(2))


.plot_diff_marg(m, type = "post_diff_ratio")
.plot_diff_marg(m, type = "post_diff_ratio", hu = TRUE)
plot_hurdle_probs(m)

