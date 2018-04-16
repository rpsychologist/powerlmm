


# Random slope SD ---------------------------------------------------------
test_that("slope_SD", {
    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "slope_SD"))

    expect_equal(get_slope_diff(p), 0.5 * 2.2 * 10)


    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_cluster_slope = 0.5,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "slope_SD"))

    expect_equal(get_slope_diff(p), 0.5 * sqrt(2.2^2 + 0.5^2) * 10)



    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "slope_SD",
                                     treatment = "treatment"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(2.2^2 + 0.5^2) * 10)

    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "slope_SD",
                                     treatment = "control"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(2.2^2) * 10)
})


# pretest_SD --------------------------------------------------------------
test_that("pretest_SD", {
    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "pretest_SD"))

    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2))


    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_cluster_slope = 0.5,
                          sigma_cluster_intercept = 2,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "pretest_SD"))

    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2^2))


    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "pretest_SD",
                                     treatment = "treatment"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2^2))

    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "pretest_SD",
                                     treatment = "control"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2))
})


# posttest_SD --------------------------------------------------------------


test_that("posttest_SD", {
    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "posttest_SD"))

    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2.2^2*10^2))


    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_cluster_slope = 0.5,
                          sigma_cluster_intercept = 2,
                          sigma_error = 10,
                          effect_size = cohend(0.5, standardizer = "posttest_SD"))

    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2.2^2*10^2 + 2^2 + 0.5^2*10^2))


    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "posttest_SD",
                                     treatment = "treatment"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2.2^2*10^2 + 2^2 + 0.5^2*10^2))

    p <- update(p, partially_nested = TRUE,
                effect_size = cohend(0.5,
                                     standardizer = "posttest_SD",
                                     treatment = "control"))
    expect_equal(get_slope_diff(p), 0.5 * sqrt(10^2 + 10^2 + 2.2^2*10^2))
})

# Raw
test_that("raw ES", {
    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = 5)

    expect_equal(get_slope_diff(p), 5)

    p <- update(p, effect_size = c(1,2,3))
    expect_equal(get_slope_diff(p), c(1,2,3))
})

# Multi Cohen's d
test_that("multi cohen's d", {
    p <- study_parameters(n1 = 11,
                          n2 = 30,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = cohend(c(0, 0.1, 0.5), standardizer = "pretest_SD"))

    expect_equal(get_slope_diff(p), c(0, 0.1, 0.5) * sqrt(10^2 + 10^2))


    p <- study_parameters(n1 = 11,
                          n2 = c(5, 2, unequal_clusters(4,5,5)),
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = cohend(c(0.1), standardizer = "pretest_SD"))

    expect_equal(get_slope_diff(p), c(0.1, 0.1, 0.1) * sqrt(10^2 + 10^2))
})



# combine types -----------------------------------------------------------
test_that("combine ES types", {
    p <- study_parameters(n1 = 11,
                          n2 = 5,
                          T_end = 10,
                          sigma_subject_intercept = 10,
                          sigma_subject_slope = 2.2,
                          sigma_error = 10,
                          effect_size = c(-5, 9,
                                          cohend(c(0.5, 0.8), standardizer = "pretest_SD"),
                                          cohend(c(0.5, 0.8), standardizer = "posttest_SD")))

    expect_equal(get_slope_diff(p), c(-5, 9, c(0.5, 0.8) * sqrt(10^2+10^2), c(0.5, 0.8) * sqrt(10^2+10^2 + 2.2^2*10^2)))
    tmp <- get_effect_size(p)
    expect_equal(tmp$ES, c(-5, 9, 0.5, 0.8, 0.5, 0.8))
    expect_equal(tmp$standardizer, c("raw", "raw", "pretest_SD", "pretest_SD", "posttest_SD", "posttest_SD"))

})


