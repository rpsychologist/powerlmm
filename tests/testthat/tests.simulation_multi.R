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

formula <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
                "wrong" = "y ~ treatment * time + (1 + time | subject)")

res <- simulate(paras, nsim = 2, formula = formula, satterthwaite = FALSE, progress = FALSE, batch_progress = FALSE)


test_that("multi_sim", {
    tmp <- summary(res)
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
    tmp <- summary(res, type = "random", para = "subject_slope")
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
    expect_error(summary(res, type = "see"), "'type' should be either 'fixed' or 'random'")
    expect_error(summary(res, para = "sdf"), "Para should be one of:")
    expect_error(summary(res, type = "random"), "No random effect named: 'time:treatment'")
    expect_error(summary(res, type = "random", para = "test", "No random effect named: 'time:treatment'"))
})




