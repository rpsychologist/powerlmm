
# Messages
test_that("ignore n3", {
    expect_message(study_parameters(n1 = 10,
                          n2 = unequal_clusters(5, 10, 20, 30),
                          n3 = per_treatment(2, 5),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5), "'n3' per_treatment argument is ignored.")

})

## Test errors


test_that("combine icc_pre_subject and sigma_cluster_intercept", {
    expect_error(study_parameters(n1 = 4,
                     n2 = 3,
                     n3 = 4,
                     sigma_cluster_intercept = 0.5,
                     icc_pre_subject = 0.5,
                     var_ratio = 0.05,
                     icc_slope = 0.1,
                     sigma_error = 0.5), "'icc_pre_subject' and 'sigma_cluster_intercept' can't be combined")

    expect_error(study_parameters(n1 = 4,
                                  n2 = 3,
                                  n3 = 4,
                                  sigma_cluster_intercept = 0.6,
                                  icc_pre_subject = 0.5,
                                  var_ratio = 0.05,
                                  icc_slope = 0.1,
                                  sigma_error = 0.5), "'icc_pre_subject' and 'sigma_cluster_intercept' can't be combined")


    expect_error(study_parameters(n1 = 4,
                                  n2 = 3,
                                  n3 = 4,
                                  sigma_cluster_intercept = c(0.4, 0.5),
                                  icc_pre_subject = c(0.5),
                                  var_ratio = 0.05,
                                  icc_slope = 0.1,
                                  sigma_error = 0.5), "'icc_pre_subject' and 'sigma_cluster_intercept' can't be combined")


})

# trying to solve for subject_slope
test_that("can't solve subject_slope", {
    msg <- "Can't use both 'icc_slope' and 'sigma_cluster_slope'"
    expect_error(study_parameters(n1 = 11,
                     n2 = 10,
                     n3 = 6,
                     icc_pre_subject = 0.5,
                     sigma_cluster_slope = 1,
                     icc_slope = 0.1,
                     cohend = -0.8), msg)

    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  sigma_cluster_slope = c(1,2),
                                  icc_slope = c(NA, 0.1),
                                  cohend = -0.8), msg)
    msg <- "'sigma_subject_slope' is missing"
    expect_error(study_parameters(n1 = 11,
                     n2 = 10,
                     n3 = 6,
                     icc_pre_subject = 0.5,
                     sigma_cluster_slope = 1,
                     cohend = -0.8), msg)

    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  sigma_cluster_slope = c(0,1),
                                  cohend = -0.8), msg)
})

# icc_pre_cluster > icc_pre_subjects
test_that("icc_pre_cluster > icc_pre_subjects", {

    msg <- "'icc_pre_cluster' can't be larger than 'icc_pre_subject'"

    # single
    expect_error(study_parameters(n1 = 3,
                                  n2 = 30,
                                  n3 = 2,
                                  T_end = 10,
                                  icc_pre_subject = c(0.1),
                                  cor_subject = -0.7,
                                  icc_pre_cluster = c(0.2),
                                  cor_cluster = 0.7,
                                  icc_slope = 0.1,
                                  var_ratio = 0.02,
                                  dropout = 0,
                                  cohend = 2
    ), msg)


    # multi
    expect_error(study_parameters(n1 = 3,
                                   n2 = 30,
                                   n3 = 2,
                                   T_end = 10,
                                   icc_pre_subject = c(0.09, 0.6),
                                   cor_subject = -0.7,
                                   icc_pre_cluster = c(0, 0.05, 0.1),
                                   cor_cluster = 0.7,
                                   icc_slope = 0.1,
                                   var_ratio = 0.02,
                                   dropout = 0,
                                   cohend = 2
    ), msg)

    # multi NA
    expect_error(study_parameters(n1 = 3,
                                  n2 = 30,
                                  n3 = 2,
                                  T_end = 10,
                                  icc_pre_subject = c(NA, 0.09, 0.6),
                                  cor_subject = -0.7,
                                  icc_pre_cluster = c(NA, 0, 0.05, 0.1),
                                  cor_cluster = 0.7,
                                  icc_slope = 0.1,
                                  var_ratio = 0.02,
                                  dropout = 0,
                                  cohend = 2
    ), msg)

})

# duplicate intercepts
test_that("cluster intercept duplicate", {
    msg <- "Can't use both 'icc_pre_cluster' and 'sigma_cluster_intercept'"
    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  n3 = 5,
                                  sigma_subject_intercept = 0.66,
                                  icc_pre_cluster = 0.25,
                                  sigma_cluster_intercept = 1.44,
                                  icc_slope = 0.05,
                                  var_ratio = 0.03,
                                  sigma_error = 1.33,
                                  cohend = 0.5), msg)

})
test_that("subject intercept duplicate", {
    msg <- "Can't use both 'icc_pre_subject' and 'sigma_subject_intercept'"
    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  n3 = 5,
                                  sigma_subject_intercept = 0.66,
                                  icc_pre_subject = 0.25,
                                  sigma_cluster_intercept = 1.44,
                                  icc_slope = 0.05,
                                  var_ratio = 0.03,
                                  sigma_error = 1.33,
                                  cohend = 0.5), msg)

})


test_that("NA and zero", {

    msg0 <- "'icc_slope' can't be > 0 when 'var_ratio' or 'sigma_subject_slope' is 0"
    msgNA <- "'icc_slope' can't be >= 0 when 'var_ratio' or 'sigma_subject_slope' is NA"
    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0,
                                  var_ratio = NA,
                                  cohend = -0.8), msgNA)

    # var_ratio NULL
    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0.05,
                                  cohend = -0.8), msgNA)

    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0.1,
                                  var_ratio = 0,
                                  cohend = -0.8), msg0)


    # multi -------------------------------------------------------------------
    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = c(0,0.5),
                                  var_ratio = c(NA, 0),
                                  cohend = -0.8), msgNA)

    # var_ratio NULL
    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = c(0, 0.05),
                                  cohend = -0.8), msgNA)

    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = c(0, 0.1),
                                  var_ratio = c(0, 1),
                                  cohend = -0.8), msg0)



   # sigma_subject_slope -----------------------------------------------------
    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0,
                                  sigma_subject_slope = NA,
                                  cohend = -0.8), msgNA)

    expect_error(study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0,
                                  sigma_subject_slope = c(NA,0,1),
                                  cohend = -0.8), msgNA)




})


test_that("duplicate information", {
    msg <- "Can't use both 'sigma_subject_slope' and 'var_ratio'"
    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  icc_pre_subject = 0.5,
                                  icc_pre_cluster = 0,
                                  sigma_subject_slope = 0.5,
                                  icc_slope = 0.05,
                                  var_ratio = 0.03,
                                  sigma_error = 1.2,
                                  cohend = 0.5), msg)

    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  icc_pre_subject = 0.5,
                                  icc_pre_cluster = 0,
                                  sigma_subject_slope = c(0, 0.5),
                                  icc_slope = 0.05,
                                  var_ratio = c(1, 0.03),
                                  sigma_error = 1.2,
                                  cohend = 0.5), msg)

})
test_that("duplicate information", {
    msg <- "Can't use 'icc_slope' with both 'sigma_subject_slope' and 'sigma_cluster_slope'"
    expect_error(study_parameters(n1 = 10,
                                   n2 = 10,
                                   n3 = 5,
                                   sigma_subject_intercept = 1.2,
                                   icc_pre_cluster = 0.05,
                                   sigma_subject_slope = 1.2,
                                   sigma_cluster_slope = 1,
                                   icc_slope = 0.05,
                                   cohend = 0.5), msg)

    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  n3 = 5,
                                  sigma_subject_intercept = 1.2,
                                  icc_pre_cluster = 0.05,
                                  sigma_subject_slope = c(0, 1.2),
                                  sigma_cluster_slope = 1,
                                  icc_slope = c(NA, 0.05),
                                  cohend = 0.5), msg)

})





# Test setup

# with only level 2 random slope
test_that("setup 3 lvl", {

    # icc_slope NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.122,
                          icc_pre_cluster = 0.1,
                          var_ratio = 0.05)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), as.numeric(NA))

    # same but icc_slope = 0
    p <- study_parameters(n1 = 10,
                        n2 = 10,
                        icc_pre_subject  = 0.122,
                        icc_pre_cluster = 0.1,
                        icc_slope = 0,
                        var_ratio =0.05)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), 0)

    # Same but with sigma_error & ICC_slope NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          var_ratio = 0.05,
                          sigma_error = 1.33)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), as.numeric(NA))
    expect_equal(p$sigma_error, 1.33)

    # Same but with sigma_error & ICC_slope 0
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          icc_slope = 0,
                          var_ratio  =0.05,
                          sigma_error = 1.33)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), 0)
    expect_equal(p$sigma_error, 1.33)
})

# With no random slopes only ICC_pre
test_that("setup 3 lvl", {

    # icc_slope & var_ratio NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject = 0.122,
                          icc_pre_cluster = 0.1)

    expect_equal(get_var_ratio(p), as.numeric(NA))
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(p$sigma_cluster_slope, NA)
    expect_equal(p$sigma_subject_slope, as.numeric(NA))


    # icc_slope NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          var_ratio = 0)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(p$sigma_cluster_slope, NA)
    expect_equal(p$sigma_subject_slope, 0)

    # same but icc_slope = 0
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          icc_slope = 0,
                          var_ratio = 0)

    expect_equal(get_var_ratio(p), 0)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(p$sigma_cluster_slope, 0)
    expect_equal(p$sigma_subject_slope, 0)

    # Same but with sigma_error & ICC_slope NULL
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          var_ratio = 0.05,
                          sigma_error = 1.33)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), as.numeric(NA))
    expect_equal(p$sigma_error, 1.33)

    # Same but with sigma_error & ICC_slope 0
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          icc_pre_subject  = 0.122,
                          icc_pre_cluster = 0.1,
                          icc_slope = 0,
                          var_ratio =0.05,
                          sigma_error = 1.33)

    expect_equal(get_var_ratio(p), 0.05)
    expect_equal(get_ICC_pre_subjects(p), 0.122)
    expect_equal(get_ICC_pre_clusters(p), 0.1)
    expect_equal(get_ICC_slope(p), 0)
    expect_equal(p$sigma_error, 1.33)


})



test_that("study setup lvl 3 minimal", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0.2,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0.2, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_cluster_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})
test_that("study setup lvl 3 minimal #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_cluster_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})
# icc_pre_cluster = NULL
test_that("study setup lvl 3 minimal icc NULL #2", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})


test_that("study setup lvl 3 minimal #3", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          icc_pre_subject = 0.2,
                          icc_pre_cluster = 0.2,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.2, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept, 0)
    expect_equal(get_ICC_pre_clusters(p), 0.2, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_cluster_intercept^2 +
                     p$sigma_error^2, 125,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})
test_that("study setup lvl 3 minimal raw", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          sigma_subject_intercept = sqrt(0.7071068^2/2),
                          sigma_subject_slope = 0.09246621,
                          sigma_cluster_intercept = sqrt(0.7071068^2/2),
                          sigma_cluster_slope = 0.0212132,
                          sigma_error = 0.7071068,
                          cohend = 0.5)
    expect_equal(get_var_ratio(p), 0.018, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0.25, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_cluster_intercept^2 +
                     p$sigma_error^2, 1,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})

test_that("study setup lvl 3 icc + sigma_error", {
    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0.2,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          sigma_error = 1.33,
                          cohend = 0.5)
    expect_equal(p$sigma_error, 1.33, tolerance = 0.001)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0.2, tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})
test_that("solve subject_intercept", {
    expect_error(study_parameters(n1 = 10,
                                  n2 = 10,
                                  n3 = 5,
                                  icc_pre_subject = 0.5,
                                  sigma_cluster_intercept = 1.44,
                                  sigma_subject_slope = 1.66,
                                  icc_slope = 0.05,
                                  var_ratio = 0.03,
                                  cohend = 0.5),
                 msg = "'icc_pre_subject' and 'sigma_cluster_intercept' can't be combined")
    # expect_equal(p$sigma_subject_slope, 1.66, tolerance = 0.001)
    # expect_equal(p$sigma_cluster_intercept, 1.44, tolerance = 0.001)
    # x <- sqrt((0.5 * (1.44^2 + (1.66^2/0.95) / 0.03))/(1 - 0.5))
    # expect_equal(p$sigma_subject_intercept, x, tolerance = 0.001)
    # expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    # expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    # expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    # icc <- 1.44^2/((1.44^2 +(1.66^2/0.95) / 0.03) * 1/0.5)
    # expect_equal(get_ICC_pre_clusters(p), icc, tolerance = 0.001)
    # expect_is(p, "plcp")
    # expect_is(p, "plcp_3lvl")
})

test_that("solve cluster_intercept", {

    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          sigma_subject_intercept = 0.66,
                          icc_pre_cluster = 0.25,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          sigma_error = 1.33,
                          cohend = 0.5)
    expect_equal(p$sigma_error, 1.33, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept, 0.66, tolerance = 0.001)
    x <- sqrt((0.25 * (1.33^2 + 0.66^2))/(1 - 0.25))
    expect_equal(p$sigma_cluster_intercept, x, tolerance = 0.001)
    expect_equal(get_var_ratio(p), 0.03, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    icc <- (0.66^2 + 0.8572242^2)/((0.66^2 + 1.33^2) * 1/0.75)
    expect_equal(get_ICC_pre_subjects(p), icc, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0.25, tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})

# common use case
test_that("solve cluster_intercept + slope", {

    p <- study_parameters(n1 = 10,
                          n2 = 10,
                          n3 = 5,
                          sigma_subject_intercept = 1.2,
                          icc_pre_cluster = 0.05,
                          sigma_subject_slope = 1.2,
                          icc_slope = 0.05,
                          sigma_error = 1.2,
                          cohend = 0.5)
    expect_equal(p$sigma_error, 1.2, tolerance = 0.001)
    expect_equal(p$sigma_subject_intercept, 1.2, tolerance = 0.001)
    x <- sqrt(0.05 * (1.2^2 + 1.2^2)/(1-0.05))
    expect_equal(p$sigma_cluster_intercept, x, tolerance = 0.001)
    x <- 1.2^2/(1-0.05) / 1.2^2
    expect_equal(get_var_ratio(p), x, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0.05, tolerance = 0.001)
    icc <- (1.2^2+0.3893314^2)/((1.2^2 + 1.2^2) * 1/0.95)
    expect_equal(get_ICC_pre_subjects(p), icc, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0.05, tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")
})


# Partially nested

test_that("Partially nested prepare", {
    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = 0.5)

    tmp <- prepare_paras(p)
    expect_equal(get_ICC_slope(tmp$control), 0)
    expect_equal(get_ICC_slope(tmp$treatment), 0.05)
})
test_that("Partially nested prepare", {
    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = 0.05,
                          var_ratio = 0.03,
                          partially_nested = TRUE,
                          cohend = 0.5)

    tmp <- prepare_paras(p)
    expect_equal(get_ICC_slope(tmp$control), 0)
    expect_equal(get_ICC_slope(tmp$treatment), 0.05)
    expect_equal(get_var_ratio(tmp$treatment), 0.03)
    expect_lt(get_var_ratio(tmp$control), 0.03)
    expect_equal(get_ICC_pre_subjects(tmp$treatment), 0.5)
    expect_equal(get_ICC_pre_subjects(tmp$control), 0.5)
    expect_equal(get_ICC_pre_clusters(tmp$treatment), 0)
    expect_equal(get_ICC_pre_clusters(tmp$control), 0)
})



test_that("multi 3lvl ICC_cluster include 0", {
    p <- study_parameters(n1 = 3,
                           n2 = 30,
                           n3 = 2,
                           T_end = 10,
                           icc_pre_subject = 0.5,
                           cor_subject = -0.7,
                           icc_pre_cluster = c(0, 0.05, 0.1),
                           cor_cluster = 0.7,
                           icc_slope = c(0, 0.1),
                           var_ratio = 0.02,
                           sigma_error = 5.5,
                           dropout = 0
    )
    expect_equal(get_ICC_slope(p), c(0,0,0, 0.1,0.1,0.1))
    expect_equal(get_var_ratio(p), rep(0.02, 6))
    expect_equal(get_ICC_pre_subjects(p), rep(0.5, 6))
    expect_equal(get_ICC_pre_clusters(p),  c(0, 0.05, 0.1, 0, 0.05, 0.1))
    expect_equal(p$sigma_error,  rep(5.5, 6))
})

test_that("multi 3lvl ICC_cluster include 0", {
    p <- study_parameters(n1 = 3,
                          n2 = 30,
                          n3 = 2,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.7,
                          icc_pre_cluster = c(0, 0.05, 0.1),
                          cor_cluster = 0.7,
                          icc_slope = c(0, 0.1),
                          var_ratio = 0.02,
                          sigma_error = 5.5,
                          dropout = 0
    )
    expect_equal(get_ICC_slope(p), c(0,0,0, 0.1,0.1,0.1))
    expect_equal(get_var_ratio(p), rep(0.02, 6))
    expect_equal(get_ICC_pre_subjects(p), rep(0.5, 6))
    expect_equal(get_ICC_pre_clusters(p),  c(0, 0.05, 0.1, 0, 0.05, 0.1))
    expect_equal(p$sigma_error,  rep(5.5, 6))
})

test_that("multi 3lvl ICC_cluster include 0, sigma_error NULL", {
    p <- study_parameters(n1 = 3,
                          n2 = 30,
                          n3 = 2,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.7,
                          icc_pre_cluster = c(0, 0.05, 0.1),
                          cor_cluster = 0.7,
                          icc_slope = c(0, 0.1),
                          var_ratio = 0.02,
                          dropout = 0
    )
    expect_equal(get_ICC_slope(p), c(0,0,0, 0.1,0.1,0.1))
    expect_equal(get_var_ratio(p), rep(0.02, 6))
    expect_equal(get_ICC_pre_subjects(p), rep(0.5, 6))
    expect_equal(get_ICC_pre_clusters(p),  c(0, 0.05, 0.1, 0, 0.05, 0.1))
    expect_equal(p$sigma_error^2, rep(100, 6))
})

# detect bug that permutations to be wrong when icc_pre_cluster included 0
test_that("multi 3lvl ICC_cluster include 0, correct grid", {
    p <- study_parameters(n1 = 3,
                          n2 = 30,
                          n3 = 2,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          cor_subject = c(0, -0.7),
                          icc_pre_cluster = c(0, 0.05),
                          cor_cluster = c(0, 0.7),
                          icc_slope = c(0.1),
                          var_ratio = 0.02,
                          dropout = 0
    )
    expect_equal(get_ICC_slope(p), rep(0.1, 8))
    expect_equal(get_var_ratio(p), rep(0.02, 8))
    expect_equal(get_ICC_pre_subjects(p), rep(0.5, 8))
    expect_equal(get_ICC_pre_clusters(p),  c(rep(0, 4), rep(0.05, 4)))
    expect_equal(p$cor_cluster, c(0, 0, 0.7, 0.7, 0, 0, 0.7, 0.7))
    expect_equal(p$cor_subject, c(0, -0.7, 0, -0.7, 0, -0.7, 0, -0.7))

})

# Test NA

test_that("study setup lvl 3 NA", {
    # no error
    p <- study_parameters(n1 = 11,
                                  n2 = 10,
                                  n3 = 6,
                                  icc_pre_subject = 0.5,
                                  icc_slope = 0,
                                  var_ratio = 0,
                                  cohend = -0.8)

    expect_equal(get_var_ratio(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")


    # no error


    p <-     study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              icc_pre_subject = 0.5,
                              icc_slope = 0.1,
                              sigma_subject_slope = 0,
                              cohend = -0.8)

    expect_equal(get_var_ratio(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 0, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")



    # no error
    p <- study_parameters(n1 = 11,
                     n2 = 10,
                     n3 = 6,
                     icc_pre_subject = 0.5,
                     sigma_cluster_slope = 1,
                     sigma_subject_slope = NA,
                     cohend = -0.8)

    expect_equal(get_var_ratio(p), 0.01, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 1, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_slope, NA)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")


    # no error
    p <- study_parameters(n1 = 11,
                     n2 = 10,
                     n3 = 6,
                     icc_pre_subject = 0.5,
                     sigma_cluster_slope = 1,
                     sigma_subject_slope = 0,
                     cohend = -0.8)

    expect_equal(get_var_ratio(p), 0.01, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), 1, tolerance = 0.001)
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_slope, 0)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")

    #
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 6,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = NA,
                          var_ratio = 0.01,
                          cohend = -0.8)

    expect_equal(get_var_ratio(p), 0.01, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), as.numeric(NA))
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")

    #
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 6,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          icc_slope = NA,
                          var_ratio = NA,
                          cohend = -0.8)

    expect_equal(get_var_ratio(p), as.numeric(NA), tolerance = 0.001)
    expect_equal(get_ICC_slope(p), as.numeric(NA))
    expect_equal(get_ICC_pre_subjects(p), 0.5, tolerance = 0.001)
    expect_equal(get_ICC_pre_clusters(p), 0)
    expect_equal(p$sigma_subject_intercept^2 +
                     p$sigma_error^2, 200,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_3lvl")


    # no error
    p <- study_parameters(n1 = 11,
                     n2 = 10,
                     n3 = 6,
                     icc_pre_subject = NA,
                     sigma_subject_slope = 1,
                     cohend = -0.8)

    expect_equal(get_var_ratio(p), 0.01, tolerance = 0.001)
    expect_equal(get_ICC_slope(p), as.numeric(NA))
    expect_equal(get_ICC_pre_subjects(p), as.numeric(NA))
    expect_equal(get_ICC_pre_clusters(p), as.numeric(NA))
    expect_equal(p$sigma_subject_slope, 1)
    expect_equal(p$sigma_cluster_slope, NA)
    expect_equal(p$sigma_subject_slope, 1)
    expect_equal(p$sigma_error^2, 100,
                 tolerance = 0.001)
    expect_is(p, "plcp")
    expect_is(p, "plcp_2lvl")

})


