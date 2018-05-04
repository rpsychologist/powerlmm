# Use for testing
set.seed(1)
n2 <- c(3, per_treatment(unequal_clusters(2, 3, 4),
                         unequal_clusters(2, 2, 3, 4)))

paras <- study_parameters(n1 = 3,
                          n2 = n2,
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.05,
                          dropout = dropout_weibull(0.3, 2),
                          cohend = -0.8)

formula <- sim_formula_compare("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
                                "wrong" = "y ~ treatment * time + (1 + time | subject)")

res <- simulate(paras, nsim = 2, formula = formula, satterthwaite = FALSE, progress = FALSE, batch_progress = FALSE)


test_that("multi_sim", {

    tmp <- summary(res, para = "treatment:time")
    expect_is(tmp, "plcp_multi_sim_summary")
    expect_output(print(tmp), "^Model: 'All' \\| Type: 'fixed'")
    expect_identical(nrow(tmp), 4L)

    tmp <- summary(res, para = "treatment:time", model = 1)
    expect_is(tmp, "plcp_multi_sim_summary")
    expect_output(print(tmp), "^Model: 'correct' \\| Type: 'fixed'")
    expect_identical(nrow(tmp), 2L)

    tmp <- summary(res, para = "treatment:time", model = "correct")
    expect_is(tmp, "plcp_multi_sim_summary")
    expect_output(print(tmp), "^Model: 'correct' \\| Type: 'fixed'")
    expect_identical(nrow(tmp), 2L)

    ## correct class for paras
    expect_true("plcp_3lvl" %in% class(res[[1]]$paras))

    ## Fixed effect

    # set 1
    x <- c(res[[1]]$res$correct$FE[4, "estimate"],
           res[[1]]$res$correct$FE[8, "estimate"])
    expect_length(x, 2)
    expect_equal(mean(x), tmp[1, "M_est"])

    # set 2
    x2 <- c(res[[2]]$res$correct$FE[4, "estimate"],
           res[[2]]$res$correct$FE[8, "estimate"])
    expect_length(x2, 2)
    expect_equal(mean(x2), tmp[2, "M_est"])

    # expect set 1 and 2 differ
    expect_true(mean(x) != mean(x2))

    ## Random effect

    # set 1
    tmp <- summary(res, para = "subject_slope", model = "correct")
    expect_identical(nrow(tmp), 2L)

    x <- c(res[[1]]$res$correct$RE[2, "vcov"],
           res[[1]]$res$correct$RE[7, "vcov"])

    expect_length(x, 2)
    expect_equal(mean(x), tmp[1, "M_est"])

    # set 2
    x2 <- c(res[[2]]$res$correct$RE[2, "vcov"],
           res[[2]]$res$correct$RE[7, "vcov"])

    expect_length(x2, 2)
    expect_equal(mean(x2), tmp[2, "M_est"])

    expect_true(mean(x2) != mean(x))

    # Test object
    expect_is(res, "plcp_multi_sim")

    ## Test saved paras
    # set 1
    x <- unlist(res[[1]]$paras$n2)
    expect_identical(x, 3)

    # set 2
    x <- res[[2]]$paras$n2
    n2 <- paras[2,"n2"]
    expect_identical(x, n2)
})


# Input validation
test_that("multi_sim summary validation", {
    # expect error
    expect_error(summary(res, para = "sdf"), "No 'para' named: sdf")
    expect_error(summary(res, model = 3), "Numeric argument 'model' is too large.")
    expect_error(summary(res, para = list("abc" = 1, "wrong" = "time"), model = 1), "'model' not found in 'para'")
    expect_error(summary(res, para = list("correct" = 1, "wrong" = "time"), model = 1), "'para' can't be a numeric value")
    expect_error(summary(res, para = list("correct" = "abc", "wrong" = "time"), model = 1), "No 'para' named: abc")
    expect_error(summary(res, para = list("correct" = "time", "wrong" = "abc")), "No 'para': abc found in 'model': wrong")
    expect_error(summary(res, para = list("correct" = "abc", "wrong" = "time", "b" = 2)), "When 'para' is a list it must be the same length as the number of models: 2")
})

test_that("sim data_transform multi", {
    p <- study_parameters(n1 = 3,
                          n2 = 5:6,
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)
    f <- sim_formula("y ~ treatment + (1 | cluster)", data_transform = transform_to_posttest, test = "treatment")
    res <- simulate(p, nsim = 2, formula = f, batch_progress = FALSE)

    x <- summary(res, para = "treatment")
    expect_is(x, "plcp_multi_sim_summary")
    expect_output(print(x), "^Model: 'All' \\| Type: 'fixed'")
    expect_output(print(x), "M_est theta")
    expect_output(print(x), "nsim:  2")

    res <- simulate(p, nsim = 2, formula = f, CI = TRUE, satterthwaite = TRUE, batch_progress = FALSE)

    x <- summary(res, para = "treatment")
    expect_output(print(x), "^Model: 'All' \\| Type: 'fixed'")
    expect_output(print(x), "M_est theta")
    expect_output(print(x), "CI_Cover CI_Wald_cover")
    expect_output(print(x), "nsim:  2")
})


