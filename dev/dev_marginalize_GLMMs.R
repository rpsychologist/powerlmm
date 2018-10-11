
# TODO
# * allow RE_level arg for post_diff plots
# * add type = "post_overlay": plot diststribution of outcomes for percentile p
# * update so plot(p_gamma, RE_level = c(2,3)) use facets instead of grid.arrange

# Gaussian ----------------------------------------------------------------
p <- study_parameters(design = study_design(),
                          n1 = 11,
                          n2 = 25,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.02,
                          effect_size = log(0.5),
                          sigma_error = 1)

m <- marginalize(p)

plot(m)
plot(m , RE = FALSE, type = "trend_dropout", RE_level = c(2))
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
plot(p_bin)
plot(m_bin , RE = FALSE, type = "trend_dropout", RE_level = c(2))
plot(m_bin, type = "post_ratio")

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
# Gamma ----------------------------------------------------------------
p_gamma <- study_parameters(design = study_design(family = "gamma"),
                         n1 = 11,
                         n2 = 25,
                         icc_pre_subject = 0.5,
                         var_ratio = 0.02,
                         effect_size = log(1),
                         shape = 1.5,
                         sigma_error = 1)

m_gamma <- marginalize(p_gamma)

plot(p_gamma, RE_level = c(2,3))
plot(m_gamma)



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

plot(m, RE = TRUE, type = "trend_dropout")
plot(m, RE = FALSE, type = "trend_dropout")

plot_hurdle_time(m$y_overall)
plot_hurdle_time(m$y_positive)
plot_hurdle_time(m$hu_prob)

.plot_diff_marg(m, type = "post_diff_ratio")
.plot_diff_marg(m, type = "post_diff_ratio", hu = TRUE)
plot_hurdle_probs(m)

