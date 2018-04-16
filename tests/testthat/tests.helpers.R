# get_icc_slope

test_that("icc_slope", {
    x <- get_ICC_slope(u1 = 1.2, v1 = 1.2)
    expect_equal(x, 0.5)

    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = 5,
                          sigma_subject_intercept = 1.55,
                          sigma_subject_slope = 1.33,
                          sigma_cluster_intercept = 1.55,
                          sigma_cluster_slope = 1.33)
    x <- get_ICC_slope(p)
    expect_equal(x, 0.5)

    # NA
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = 5,
                          sigma_subject_intercept = 1.55,
                          sigma_subject_slope = NA,
                          sigma_cluster_intercept = 1.55,
                          sigma_cluster_slope = 1)
    x <- get_ICC_slope(p)
    expect_equal(x, 1)

    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = 5,
                          sigma_subject_intercept = 1.55,
                          sigma_subject_slope = 1,
                          sigma_cluster_intercept = 1.55,
                          sigma_cluster_slope = NA)
    x <- get_ICC_slope(p)
    expect_equal(x, as.numeric(NA))

    # Multi
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = 5,
                          sigma_subject_intercept = 1.55,
                          sigma_subject_slope = c(1.33, sqrt(1.33^2*2)),
                          sigma_cluster_intercept = 1.55,
                          sigma_cluster_slope = c(1.33, sqrt(1.33^2*2)))

    x <- get_ICC_slope(p)
    expect_equal(x, c(0.5, 1/3, 2/3, 0.5))


    })

# get_var_ratio
test_that("var_ratio", {
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = 1.2,
                              sigma_subject_slope = 0.2,
                              sigma_cluster_intercept = 0,
                              sigma_cluster_slope = 0.25,
                              sigma_error = 1.2,
                              cohend = -0.8)

    x <- get_var_ratio(paras)
    expect_equal(x, (0.2^2 + 0.25^2)/(1.2^2))

    x <- get_var_ratio(u1 = 0.222, v1 = 0.333, error = 1.33)
    expect_equal(x, (0.222^2 + 0.333^2)/(1.33^2))

    # Multi
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = 1.2,
                              sigma_subject_slope = c(sqrt(0.6), sqrt(0.3)),
                              sigma_cluster_intercept = 0,
                              sigma_cluster_slope = c(sqrt(0.6), sqrt(0.3)),
                              sigma_error = c(sqrt(1.2), sqrt(2.4)),
                              cohend = -0.8)

    x1 <- get_var_ratio(paras)
    x2 <- with(paras, (sigma_subject_slope^2 + sigma_cluster_slope^2)/sigma_error^2)
    expect_equal(x1, x2)
    expect_equal(x1, c(1, 3/4, 3/4, 1/2, 1/2, 3/8, 3/8, 1/4))
})

# get_ICC_pre_subjects
test_that("ICC_pre_subjects", {
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = 1.2,
                              sigma_subject_slope = 0.2,
                              sigma_cluster_intercept = 0.5,
                              sigma_cluster_slope = 0.2,
                              sigma_error = 1.25,
                              cohend = -0.8)

    x <- get_ICC_pre_subjects(paras)
    expect_equal(x, (1.2^2 + 0.5^2)/(1.2^2 + 0.5^2 + 1.25^2))

    x <- get_ICC_pre_subjects(u0 = 1.33, v0 = 0.5, error = 1.55)
    expect_equal(x, (1.33^2+ 0.5^2)/(1.33^2 + 0.5^2 + 1.55^2))


    # Multi
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = c(1.22, 2.33, 4.55),
                              sigma_subject_slope = 0.2,
                              sigma_cluster_intercept = c(0, 0.5),
                              sigma_cluster_slope = 0.2,
                              sigma_error = c(1.25, 0.55),
                              cohend = -0.8)

    expect_equal(nrow(paras), 3*2*2)
    x1 <- get_ICC_pre_subjects(paras)
    x2 <- with(paras, (sigma_subject_intercept^2 + sigma_cluster_intercept^2)/(sigma_subject_intercept^2 +
                                                     sigma_cluster_intercept^2 +
                                                     sigma_error^2))
    expect_equal(x1, x2)
})

# get_ICC_pre_clusters
test_that("ICC_pre_clusters", {
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = 1.2,
                              sigma_subject_slope = 0.2,
                              sigma_cluster_intercept = 0.5,
                              sigma_cluster_slope = 0.2,
                              sigma_error = 1.25,
                              cohend = -0.8)

    x <- get_ICC_pre_clusters(paras)
    expect_equal(x, 0.5^2/(1.2^2 + 0.5^2 + 1.25^2))

    x <- get_ICC_pre_clusters(u0 = 1.33, v0 = 0.5, error = 1.55)
    expect_equal(x, 0.5^2/(1.33^2 + 0.5^2 + 1.55^2))


    # Multi
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 3,
                              T_end = 10,
                              sigma_subject_intercept = c(1.22, 2.33),
                              sigma_subject_slope = 0.2,
                              sigma_cluster_intercept = c(0, 0.5, 1.2),
                              sigma_cluster_slope = 0.2,
                              sigma_error = c(1.25, 0.55),
                              cohend = -0.8)

    expect_equal(nrow(paras), 3*2*2)
    x1 <- get_ICC_pre_clusters(paras)
    x2 <- with(paras, sigma_cluster_intercept^2/(sigma_subject_intercept^2 +
                                                     sigma_cluster_intercept^2 +
                                                     sigma_error^2))
    expect_equal(x1, x2)
})



## Testers

# is.unequal_clusters and per_treatment
test_that("is.unequal_clusters per_treatment", {
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(1,2,3,4),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)

    expect_true(is.unequal_clusters(p$n2))

    # unequal in 1 arm
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(1,2,3,4), 3),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)
    expect_true(is.unequal_clusters(p$n2))
    expect_true(is.per_treatment(p$n2))

    # different unequal in each arms
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(1,2,3,4),
                                             unequal_clusters(10,15,3,4)),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)
    expect_true(is.unequal_clusters(p$n2))
    expect_true(is.per_treatment(p$n2))

    # no unequal clusters
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)
    expect_false(is.unequal_clusters(p$n2))
    expect_false(is.per_treatment(p$n2))

    # no unequal, both n2 per treatment
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(10, 15),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)
    expect_false(is.unequal_clusters(p$n2))
    expect_true(is.per_treatment(p$n2))

    # per_treatment n2 n3
    # no unequal, both n2 per treatment
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(10, 15),
                          n3 = per_treatment(5,8),
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.05,
                          cohend = -0.8)
    expect_false(is.unequal_clusters(p$n2))
    expect_true(is.per_treatment(p$n2))

})



# SDS ---------------------------------------------------------------------
test_that("get_SDS", {
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          sigma_subject_intercept = 1L,
                          icc_slope = 0,
                          sigma_error = 1L,
                          icc_pre_cluster = 0L,
                          var_ratio = 0L)

    x <- get_sds(p)
    expect_identical(nrow(x), 11L)
    expect_equal(x$SD_with_random_slopes, rep(sqrt(1^2 + 1^2), 11))
    expect_error(x, NA)

    # varying slope
    p <- update(p, var_ratio = 0.05, icc_slope = 0.1)
    x <- get_sds(p)
    expect_identical(nrow(x), 11L)
    expect_equal(x$SD_no_random_slopes, rep(sqrt(1^2 + 1^2), 11))

    expect_length(unique(x$SD_with_random_slopes), 11)

    expect_error(plot(x), NA)

    # NA
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          sigma_subject_intercept = 1L,
                          icc_slope = 0,
                          sigma_error = 1L,
                          icc_pre_cluster = NA,
                          var_ratio = 0L)

    x <- get_sds(p)
    expect_identical(nrow(x), 11L)
    expect_equal(x$SD_with_random_slopes, rep(sqrt(1^2 + 1^2), 11))
    expect_error(x, NA)


    # Control, partially nested
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(20, 40),
                          n3 = 4,
                          sigma_subject_intercept = 10,
                          sigma_cluster_intercept = 5,
                          sigma_error = 10,
                          partially_nested = TRUE)

    expect_true(all(get_sds(p, treatment = "treatment")$SD_with_random_slopes == sqrt(10^2 + 5^2 + 10^2)))
    expect_true(all(get_sds(p, treatment = "control")$SD_with_random_slopes == sqrt(10^2 + 10^2)))

})

# get_correlation_matrix --------------------------------------------------
test_that("get_correlation_matrix", {
    p1 <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          sigma_subject_intercept = 1L,
                          icc_slope = 0,
                          sigma_error = 1L,
                          icc_pre_cluster = 0L,
                          var_ratio = 0L)

    x <- get_correlation_matrix(p1)
    tmp <- dplyr::near(x[lower.tri(x)], 0.5)

    expect_true(all(tmp))
    expect_true(all(diag(x) == 1))
    expect_error(plot(x), NA)

    # NA
    p1 <- study_parameters(n1 = 11,
                           n2 = 10,
                           n3 = 3,
                           T_end = 10,
                           sigma_subject_intercept = 1L,
                           icc_slope = 0,
                           sigma_error = 1L,
                           icc_pre_cluster = NA,
                           var_ratio = 0L)

    x <- get_correlation_matrix(p1)
    tmp <- dplyr::near(x[lower.tri(x)], 0.5)

    expect_true(all(tmp))
    expect_true(all(diag(x) == 1))
    expect_error(plot(x), NA)

    # varying slope
    p1 <- update(p1, var_ratio = 0.05, icc_slope = 0.1)
    x <- get_correlation_matrix(p1)
    expect_true(all(diag(x) == 1))
    expect_equal(dim(x), c(11,11))
    expect_true(all(!is.na(x[lower.tri(x)])))

})
# get_VPC --------------------------------------------------
test_that("get_VPC", {
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          sigma_subject_intercept = 1L,
                          icc_slope = 0,
                          sigma_error = 1L,
                          icc_pre_cluster = 0L,
                          var_ratio = 0L)

    x <- get_VPC(p)
    expect_identical(nrow(x), 11L)
    expect_equal(x$between_clusters, rep(0, 11))
    expect_equal(x$between_subjects, rep(50, 11))
    expect_equal(x$within_subjects, rep(50, 11))
    expect_equal(x$tot_var, rep(0, 11))
    expect_error(x, NA)

    # NA
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          sigma_subject_intercept = 1L,
                          icc_slope = 0,
                          sigma_error = 1L,
                          icc_pre_cluster = NA,
                          var_ratio = 0L)

    x <- get_VPC(p)
    expect_identical(nrow(x), 11L)
    expect_equal(x$between_clusters, rep(0, 11))
    expect_equal(x$between_subjects, rep(50, 11))
    expect_equal(x$within_subjects, rep(50, 11))
    expect_equal(x$tot_var, rep(0, 11))
    expect_error(x, NA)

    # varying slope
    p <- update(p, var_ratio = 0.05, icc_slope = 0.1)
    x <- get_VPC(p)
    expect_identical(nrow(x), 11L)
    expect_length(unique(x$between_clusters), 11)
    expect_error(plot(x), NA)

})
