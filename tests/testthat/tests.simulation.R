p <- study_parameters(n1 = 10,
                      n2 = 10,
                      n3 = 5,
                      sigma_subject_intercept = 1.44,
                      icc_pre_cluster = 0,
                      sigma_subject_slope = 0.2,
                      icc_slope = 0.05,
                      sigma_error = 1.44,
                      cohend = 0.5)


# munge_results -----------------------------------------------------------


test_that("munge_results", {
    set.seed(5)

    formula <- list("correct" = "y ~ treatment * time + (1 + time | subject) +
                    (0 + time | cluster)")
    res <- lapply(1:3, simulate_,
                  paras = p,
                  satterthwaite = FALSE,
                  CI = FALSE,
                  formula = formula)
    munged <- munge_results(list(res = res))

    # subject_slope
    tmp <- rep(NA, 3)
    for(i in 1:3) {
        x <- res[[i]]$correct$RE
        tmp[i] <- x[x$parameter == "subject_time", "vcov"]
    }
    x <- munged$res$correct$RE
    x <- x[x$parameter == "subject_slope", "vcov"]
    expect_equal(x, tmp)

    # cluster_slope
    tmp <- rep(NA, 3)
    for(i in 1:3) {
        x <- res[[i]]$correct$RE
        tmp[i] <- x[x$parameter == "cluster_time", "vcov"]
    }
    x <- munged$res$correct$RE
    x <- x[x$parameter == "cluster_slope", "vcov"]
    expect_equal(x, tmp)

    # time:treatment
    # subject_slope
    tmp <- rep(NA, 3)
    for(i in 1:3) {
        x <- res[[i]]$correct$FE
        tmp[i] <- x[x$parameter == "time:treatment", "estimate"]
    }
    x <- munged$res$correct$FE
    x <- x[x$parameter == "time:treatment", "estimate"]
    expect_equal(x, tmp)

})


# extract_results ---------------------------------------------------------
test_that("extract results", {
    set.seed(34534)
    d <- simulate_data(p)
    fit <- lme4::lmer(y ~ treatment * time + (1 + time | subject) +
                          (0 + time | cluster), data = d)
    tmp <- extract_results(list(fit), CI = FALSE, paras = p)

    x <- tmp[[1]]$RE
    expect_equal(x$vcov, c(2.330233, 0.032569, 0.008725, 2.070966, -0.193),
                 tolerance = 0.001)

    pnames <- c( "subject_(Intercept)", "subject_time", "cluster_time",
                 "Residual_NA", "subject_(Intercept)_time")
    expect_equal(x$parameter, pnames)
})
test_that("extract results satterthwaite", {
    set.seed(34534)
    d <- simulate_data(p)
    fit <- lmerTest::lmer(y ~ treatment * time + (1 + time | subject) +
                              (0 + time | cluster), data = d)
    tmp <- extract_results(list(fit), CI = FALSE, paras = p)

    x <- tmp[[1]]$RE
    expect_equal(x$vcov, c(2.330233, 0.032569, 0.008725, 2.070966, -0.193),
                 tolerance = 0.001)

    pnames <- c( "subject_(Intercept)", "subject_time", "cluster_time",
                 "Residual_NA", "subject_(Intercept)_time")
    expect_equal(x$parameter, pnames)
})



test_that("simulation summary", {

    set.seed(5446)
    formula <- list("correct" = "y ~ treatment * time + (1 + time | subject) +
                    (0 + time | cluster)")
    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- summary(res)

    # Est
    est <- c(res$res$correct$FE[c(4,8,12), "estimate"])
    y <- mean(est)
    expect_equal(tmp$summary$correct$FE[4, "M_est"], y, tolerance = 0.001)

    # SE
    se <- c(res$res$correct$FE[c(4,8,12), "se"])
    y <- mean(se)
    expect_equal(tmp$summary$correct$FE[4, "M_se"], y, tolerance = 0.001)

    # Emperical SE
    expect_equal(tmp$summary$correct$FE[4,"SD_est"], sd(est), tolerance = 0.001)


    # Subject_slope
    est <- res$res$correct$RE
    est <- est[est$parameter == "subject_slope", "vcov"]

    expect_equal(tmp$summary$correct$RE[2, "M_est"], mean(est), tolerance = 0.001)
    expect_equal(tmp$summary$correct$RE[2, "prop_zero"], mean(abs(est - 0) < 0.00001),
                 tolerance = 0.001)

    # Subject intercept_slope cor
    est <- res$res$correct$RE
    est <- est[est$parameter == "cor_subject", "vcov"]
    expect_equal(tmp$summary$correct$RE[5, "M_est"], mean(est), tolerance = 0.001)


    # Cluster slope
    est <- res$res$correct$RE
    est <- est[est$parameter == "cluster_slope", "vcov"]

    expect_equal(tmp$summary$correct$RE[3, "M_est"], mean(est), tolerance = 0.001)
    expect_equal(tmp$summary$correct$RE[3, "prop_zero"], mean(abs(est - 0) < 0.00001),
                 tolerance = 0.001)

})

test_that("simulation summary alpha", {

    set.seed(5446)
    formula <- list("correct" = "y ~ treatment * time + (1  | subject)")
    res <- simulate(p, nsim = 10, formula = formula, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- summary(res)[[1]]$correct$FE[, "Power"]

    tmp2 <- summary(res, alpha = 0.5)[[1]]$correct$FE[, "Power"]

    expect_true(all(tmp < tmp2))
})

test_that("simulation partially nested", {
    set.seed(45443)
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 4,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          cor_subject = -0.5,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = -0.5)
    ##
    formula <- "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)"

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    ##
    formula <- "y ~ treatment * time + (1 + time | subject) + (1 + treatment:time | cluster)"

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    ##
    formula <- "y ~ treatment * time + (1 + time | subject) + (0 + time:treatment | cluster)"

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    ##
    formula <- "y ~ treatment * time + (1 + time | subject) + (1 + time:treatment | cluster)"

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

})

# Test formula
test_that("Simulation formula (character)", {

    expect_is(check_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"), "list")
    expect_error(check_formula("y ~ treatment2 * time + (1 + time | subject) + (0 + time | cluster)"))
    expect_is(check_formula("y ~ treatment + time + treatment:time + (1 + time | subject) + (0 + time | cluster)"), "list")
    expect_error(check_formula("y ~ treatment + time + treatment:time + (1 + time2 | subject) + (0 + time | cluster)"))
    expect_error(check_formula("log(y) ~ treatment + time + treatment:time + (1 + time | subject) + (0 + time | cluster)"))
    expect_error(check_formula("log(y) ~ treatment + poly(time, 2) + treatment:time + (1 + time | subject) + (0 + time | cluster)"))

})
# Test formula
test_that("Simulation formula (list)", {


    expect_is(check_formula(list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")), "list")
    expect_is(check_formula(list("wrong"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")), "list")

    f <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "wrong" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_is(check_formula(f), "list")

    # same names
    f <- list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Both formulas can't have the same name")

    # Wrong names
    f <- list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "Wong"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Formula names must be either 'correct' or 'wrong'")

    # No names
    f <- list("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Formula should be a named list")

    # Wrong terms
    f <- list("correct" = "y ~ treatment2 * time + (1 + time | subject) + (0 + time | cluster)",
              "wrong" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "treatment2 is not an allowed variable name.")

})




# Test simulations run without error
test_that("Simulation runs, unequal_clusters", {
    p <- study_parameters(n1 = 3,
                          n2 = unequal_clusters(15,15),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)
    f <- "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"
    res <- simulate(p, nsim = 2, formula = f, satterthwaite = FALSE,
                    progress = FALSE)

    expect_is(res$res$correct$FE[1,2], "numeric")

})
test_that("Simulation runs, dropout", {
    p <- study_parameters(n1 = 5,
                          n2 = unequal_clusters(15,15),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)
    f <- "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"
    res <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE,
                    progress = FALSE)

    expect_is(res$res$correct$FE[1,2], "numeric")

})






