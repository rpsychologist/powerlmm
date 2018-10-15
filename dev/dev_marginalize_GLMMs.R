
# TODO
# * allow RE_level arg for post_diff plots
# * add type = "post_overlay": plot diststribution of outcomes for percentile p
# * update so plot(p_gamma, RE_level = c(2,3)) use facets instead of grid.arrange
# * add plote type "trend_dist": plot RE dist at each time point, like ridge plots
#    - should be possible to overlay tx and cc in same panel?
# * add some type of level 1 plot of resp distribution, e.g. for different combinations of RE effects?
# * allow comparing subject at different cluster level RE effects?

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
plot(m, type = "trend_ridges", RE_level = c(2,3))

# TODO: support link scale
plot(p, type = "trend_ridges", RE_level = c(2,3))

plot(m , RE = FALSE, type = "post_diff", RE_level = c(2))
plot(m , RE = FALSE, type = "post_diff", RE_level = c(3))

# Binomial ----------------------------------------------------------------
p_bin <- study_parameters(design = study_design(family = "binomial"),
                          n1 = 11,
                          n2 = 25,
                          fixed_intercept = qlogis(0.7),
                          sigma_subject_intercept = 1,
                          sigma_subject_slope = 0,
                          cor_subject = -0.5,
                          effect_size = log(0.5),
                          sigma_error = 1)

m_bin <- marginalize(p_bin, hu = TRUE)

plot(m_bin, RE_level = c(1, 2))


plot(m_bin,
     RE = FALSE,
     type = "trend",
     RE_level = c(2,3))

plot(m_bin,
     RE = FALSE,
     type = "trend_ridges",
     RE_level = c(2,3))

plot(m_bin, type = "post_ratio")
plot(m_bin, type = "post_diff")


# Poisson ----------------------------------------------------------------
p_pois <- study_parameters(design = study_design(family = "poisson"),
                          n1 = 11,
                          n2 = 25,
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
                         n1 = 11,
                         n2 = 25,
                         fixed_intercept = log(500),
                         sigma_subject_intercept = 1,
                         sigma_cluster_intercept = 0.2,
                         effect_size = log(0.5),
                         shape = 2,
                         sigma_error = 1)

m_gamma <- marginalize(p_gamma)

plot(p_gamma, RE_level = c(2,3))
plot(m_gamma)

plot(m_gamma, type = "trend", RE_level = c(1,2), sd2_p = c(0.5, 0.5)) + ylim(0, 5000) + scale_y_log10()

plot(m_gamma, type = "trend_ridges", RE_level = c(1,2)) + xlim(0, 3000) + scale_x_log10()


plot(m_gamma, type = "trend_ridges", RE_level = c(2,3)) + xlim(0, 5000)


# Hurdle models


# CTP hurdle lognormal ----------------------------------------------------
p <- study_parameters(
    design = study_design(family = "hurdle"),
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
    RR_cont = 0.33,
    OR_hu = 2,
    marginal = FALSE,
    family = "gamma")

m <- marginalize(p)

# Done
plot(m, RE = FALSE)
plot(m, RE = TRUE)

plot(m, RE = TRUE, type = "trend")
plot(m, RE = FALSE, type = "trend")

# TODO: level 1
plot(m, RE = TRUE, type = "trend_ridges", RE_level = c(2))

library(gridExtra)

# TODO: fix so this works. plot should return ggplot objects
grid.arrange(plot(m, RE = FALSE, type = "post_diff"), plot(m, RE = FALSE, type = "post_ratio"))



plot_hurdle_time(m$y_overall)
plot_hurdle_time(m$y_positive)
plot_hurdle_time(m$hu_prob)

.plot_diff_marg(m, type = "post_diff_ratio")
.plot_diff_marg(m, type = "post_diff_ratio", hu = TRUE)
plot_hurdle_probs(m)

