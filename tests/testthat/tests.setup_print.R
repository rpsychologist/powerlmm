# Some tests that print.plcp is working


# two level ---------------------------------------------------------------
test_that("setup print", {

    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = NA,
                          sigma_subject_slope = 0.2,
                          icc_slope = NA,
                          sigma_error = 1.44,
                          cohend = 0.5)

    expect_output(str(print(p)),  "Study setup \\(two-level\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = 10000")
    expect_output(str(print(p)), "10010 \\(total\\)")

    #  partially nested will do nothing
    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = NA,
                          sigma_subject_slope = 0.2,
                          icc_slope = NA,
                          sigma_error = 1.44,
                          partially_nested =TRUE,
                          cohend = 0.5)

    expect_output(str(print(p)), "Study setup \\(two-level\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = 10000")
    expect_output(str(print(p)), "10010 \\(total\\)")

    # with dropout
    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = NA,
                          sigma_subject_slope = 0.2,
                          icc_slope = NA,
                          sigma_error = 1.44,
                          partially_nested =TRUE,
                          dropout = dropout_weibull(.3, 2),
                          cohend = 0.5)

    expect_output(str(print(p)),  "Study setup \\(two-level\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = 10000")
    expect_output(str(print(p)), "10010 \\(total\\)")
    expect_output(str(print(p)), " 0,  0,  2,  4,  7, 10, 15, 19, 25, 30 \\(%, control\\)")
})




# three-level -------------------------------------------------------------
test_that("setup print", {

    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          n3 = 5,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = NA,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)

    expect_output(str(print(p)),  "Study setup \\(three-level\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = 10000 x 5")
    expect_output(str(print(p)), "n3 = 5")
    expect_output(str(print(p)), "50050")

    #  partially nesting, total_n should not chabnge
    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          n3 = 5,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = 0.5)

    expect_output(str(print(p)),  "Study setup \\(three-level, partially nested\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = 10000 x 5")
    expect_output(str(print(p)), "n3 = 5")
    expect_output(str(print(p)), "50050")

    # with dropout
    p <- study_parameters(n1 = 10,
                          n2 = per_treatment(10, 10000),
                          n3 = 5,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = NA,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(.3, 2),
                          cohend = 0.5)

    expect_output(str(print(p)),  "Study setup \\(three-level\\)")
    expect_output(str(print(p)), "n2 = 10000 x 5")
    expect_output(str(print(p)), "50050")
    expect_output(str(print(p)), " 0,  0,  2,  4,  7, 10, 15, 19, 25, 30 \\(%, control\\)")
})

# Random n2
test_that("setup print", {

    #  partially nesting, total_n should not chabnge
    p <- study_parameters(n1 = 10,
                          n2 = unequal_clusters(func = rpois(6, 5)),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = 0.5)

    expect_output(str(print(p)), "Study setup \\(three-level, partially nested\\)")
    expect_output(str(print(p)), "n1 = 10")
    expect_output(str(print(p)), "n2 = rpois\\(6, 5\\)")
    expect_output(str(print(p)), "n3 = 6")
})
