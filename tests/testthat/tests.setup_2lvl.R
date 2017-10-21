test_that("study setup lvl 2 minimal", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.03,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 + p$sigma_error^2, 1,
                 tolerance = 0.001)
    expect_is(p, "plcp")
})
test_that("study setup lvl 2 minimal #unstandardized", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_intercept = 1.2,
                          sigma_subject_slope = 0.5,
                          sigma_error = 1.2,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.1736111, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 + p$sigma_error^2, 1.2^2+1.2^2,
                 tolerance = 0.001)
    expect_is(p, "plcp")
})

# expect error
test_that("study setup lvl 2 error multiple slope variance", {
    msg <- "'sigma_subject_slope' or 'var_ratio' or 'sigma_error' should be NULL"
    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  sigma_subject_slope = 1,
                                  var_ratio = 0.03,
                                  sigma_error = 1), msg)

})
test_that("missing sigma_subject_intercept and icc_pre", {
    msg <- "Both 'sigma_subject_intercept' and 'icc_pre_subject'"
    expect_error(study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_error = 0.7071068,
                          var_ratio = 0.03), msg)
})
test_that("solve with 0 'var_ratio'", {
    msg <- "'var_ratio' can't be zero"
    expect_error(study_parameters(n1 = 10,
                                   n2 = 10,
                                   sigma_subject_slope = 0.1224745,
                                   icc_pre_subject = 0.2,
                                   var_ratio = 0), msg)
})


# Tests

#Compound symmetry 2lvl
test_that("setup 2 lvl cs", {

    # icc_slope & var_ratio NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.5)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)


    # icc_slope NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          var_ratio = 0)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)

    # same but icc_slope = 0
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_slope = 0,
                          var_ratio = 0)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)

    # same but icc_slope & icc_pre_cluster= 0
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0,
                          var_ratio = 0)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)

    ## Add sigma_error
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0,
                          var_ratio = 0,
                          sigma_error = 1.2)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_error, 1.2)

    # remove icc_pre_cluster
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_slope = 0,
                          var_ratio = 0,
                          sigma_error = 1.2)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_error, 1.2)

    # remove var_ratio
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0,
                          sigma_error = 1.2)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_error, 1.2)

    # remove icc_slope
    # add var_ratio
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0,
                          sigma_error = 1.2)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_error, 1.2)

    # remove all
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.5,
                          sigma_error = 1.2)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.5)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_error, 1.2)




})


test_that("study setup lvl 2 solve slope", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})

test_that("study setup lvl 2 solve slope #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.5,
                          sigma_error = 0.7071068,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})
test_that("study setup lvl 2 solve slope #3", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_intercept = 0.7071068,
                          sigma_error = 0.7071068,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})


# solve for intercept var
test_that("solve sigma_intercept & sigma_error", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.4,
                          sigma_subject_slope = 0.1224745,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.4, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})
test_that("solve sigma_intercept & sigma_slope", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.4,
                          sigma_error =  0.7071068,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.4, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})
test_that("solve sigma_intercept #3", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.4,
                          sigma_error =  0.7071068,
                          sigma_subject_slope = 0.1224745)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.4, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})
test_that("solve sigma_error", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          sigma_subject_intercept = 0.5773503,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.4, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})

test_that("solve sigma_error #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          icc_pre_subject = 0.4,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.4, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_is(p, "plcp")
})


# multiple
test_that("multi 'icc_pre_subject'", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          icc_pre_subject = c(0, 0.1, 0.4),
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), c(0.03, 0.03, 0.03), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0, 0.1, 0.4), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0), tolerance = 0.001)
    expect_identical(nrow(p), 3L)
    expect_is(p, "plcp_multi")
})
test_that("multi 'var_ratio'", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          icc_pre_subject = 0.22,
                          var_ratio = c(0.01, 0.03, 0.51))
    expect_equal(get_var_ratio(p), c(0.01, 0.03, 0.51), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.22, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0), tolerance = 0.001)
    expect_identical(nrow(p), 3L)
    expect_is(p, "plcp_multi")
})
test_that("multi 'var_ratio' with 0", {
    p <- study_parameters(n1 = 11,
                           n2 = 30,
                           icc_pre_subject = 0.5,
                           var_ratio = c(0, 0.1, 0.2),
                           cohend = -0.5)
    expect_equal(get_var_ratio(p), c(0, 0.1, 0.2))
    expect_equal(get_ICC_slope(p), c(NA, 0, 0))
    expect_equal(get_ICC_pre_subjects(p), c(0.5, 0.5, 0.5))

})

test_that("multi 'sigma_subject_slope'", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = c(0.1224745, 0.3224745, 4),
                          icc_pre_subject = 0.22,
                          var_ratio = 0.03)
    expect_equal(get_var_ratio(p), c(0.03, 0.03, 0.03), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.22, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0), tolerance = 0.001)
    expect_identical(nrow(p), 3L)
    expect_is(p, "plcp_multi")
})

test_that("multi 'sigma_error'", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          icc_pre_subject = 0.22,
                          sigma_error = c(0.5, 1.2, 2.3))
    expect_equal(get_var_ratio(p), c(0.06000001,
                                     0.01041667,
                                     0.002835539), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.22, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0), tolerance = 0.001)
    expect_identical(nrow(p), 3L)
    expect_is(p, "plcp_multi")
})
test_that("multi 'sigma_error' #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = 0.1224745,
                          sigma_subject_intercept = c(0.3, 1.2, 3.3),
                          var_ratio = 0.033)
    expect_equal(get_var_ratio(p), c(0.033,
                                     0.033,
                                     0.033), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.1652754,
                                            0.7600767,
                                            0.9599327), tolerance = 0.001)
    expect_equal(p$sigma_error, c(0.6741999,
                                  0.6741999,
                                  0.6741999), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0), tolerance = 0.001)
    expect_identical(nrow(p), 3L)
    expect_is(p, "plcp_multi")
})

# combine multi
test_that("multi 'sigma_subject_slope'", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          sigma_subject_slope = c(0.1224745),
                          icc_pre_subject = c(0.1, 0.22),
                          var_ratio = c(0.01, 0.03))
    expect_equal(get_var_ratio(p), c(0.01, 0.03, 0.01, 0.03), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.1, 0.1, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(p$sigma_subject_slope, rep(0.1224745, 4), tolerance = 0.001)
    expect_identical(nrow(p), 4L)
    expect_is(p, "plcp_multi")
})

# combine multi
test_that("multi with unequal_clusters", {
    p <- study_parameters(n1 = 10,
                          n2 = unequal_clusters(4, 10, 25),
                          sigma_subject_slope = c(0.1224745),
                          icc_pre_subject = c(0.1, 0.22),
                          var_ratio = c(0.01, 0.03))
    expect_equal(sum(unlist(p$n2)), (4+10+25)*4)
    expect_equal(get_var_ratio(p), c(0.01, 0.03, 0.01, 0.03), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.1, 0.1, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(p$sigma_subject_slope, rep(0.1224745, 4), tolerance = 0.001)
    expect_identical(nrow(p), 4L)
    expect_is(p, "plcp_multi")
})
test_that("multi with dropout #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 5,
                          sigma_subject_slope = c(0.1224745),
                          icc_pre_subject = c(0.1, 0.22),
                          var_ratio = c(0.01, 0.03),
                          dropout = dropout_weibull(0.1, 2))
    expect_equal(get_ICC_slope(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.1, 0.1, 0.22, 0.22), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(p$sigma_subject_slope, rep(0.1224745, 4), tolerance = 0.001)
    expect_identical(nrow(p), 4L)
    expect_is(p, "plcp_multi")
})

test_that("multi dropout", {
    p <- study_parameters(n1 = 10,
                          n2 = 5,
                          sigma_subject_slope = c(0.1224745),
                          icc_pre_subject = 0.4,
                          var_ratio = c(0.03, 0.1),
                          dropout = c(dropout_weibull(0.1, 2),
                                      dropout_weibull(0.4, 2)))
    expect_equal(get_ICC_slope(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), c(0.4, 0.4, 0.4, 0.4), tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), c(0, 0, 0, 0), tolerance = 0.001)
    expect_equal(p$sigma_subject_slope, rep(0.1224745, 4), tolerance = 0.001)
    expect_equal(p$dropout[[1]](0:9)[10], 1-0.1)
    expect_equal(p$dropout[[2]](0:9)[10], 1-0.1)
    expect_equal(p$dropout[[3]](0:9)[10], 1-0.4)
    expect_equal(p$dropout[[4]](0:9)[10], 1-0.4)
    expect_identical(nrow(p), 4L)
    expect_is(p, "plcp_multi")
})

# Multi n1 with T_end default
test_that("multi n1 no T_end", {
    p <- study_parameters(n1 = 5:6,
                          n2 = 6:7,
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)

    expect_equal(nrow(p), 4)
    expect_equal(p$T_end, c(4,5,4,5))
    expect_equal(p$n1, c(5,6,5,6))
})


