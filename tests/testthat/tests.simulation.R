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

    formula <- compare_sim_formulas("correct" = sim_formula("y ~ treatment * time + (1 + time | subject) +
                    (0 + time | cluster)"))
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
        tmp[i] <- x[x$parameter == "treatment:time", "estimate"]
    }
    x <- munged$res$correct$FE
    x <- x[x$parameter == "treatment:time", "estimate"]
    expect_equal(x, tmp)

})


# extract_results ---------------------------------------------------------
test_that("extract results", {
    set.seed(34534)
    d <- simulate_data(p)
    fit <- lme4::lmer(y ~ treatment * time + (1 + time | subject) +
                          (0 + time | cluster), data = d)
    tmp <- extract_results(list(list("fit" = fit)), CI = FALSE, df_bw = 8, tot_n = 100, sim = 1)

    x <- tmp[[1]]$RE
    expect_equal(x$vcov, c(2.330233, 0.032569, 0.008725, 2.070966, -0.193),
                 tolerance = 0.001)

    pnames <- c( "subject_(Intercept)", "subject_time", "cluster_time",
                 "Residual_NA", "subject_(Intercept)_time")
    expect_equal(x$parameter, pnames)

    # satterth
    tmp <- extract_results(list(list("fit" = fit, "test" = "treatment:time")), satterthwaite = TRUE, CI = FALSE, df_bw = 8, tot_n = 100, sim = 1)

    expect_false(is.na(tmp[[1]]$FE[4,"pval"]))
    expect_false(is.na(tmp[[1]]$FE[4,"df"]))


})


test_that("simulation summary", {

    p <- update(p, fixed_intercept = 4.4,
                   fixed_slope = -0.22)

    set.seed(5446)
    formula <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- summary(res)

    # params
    x <- as.character(tmp$summary$default$FE$parameter)
    expect_equal(x, c("(Intercept)", "treatment", "time", "treatment:time"))

    # theta
    expect_equal(tmp$summary$default$FE$theta, c(4.4, 0, -0.22, 0.1131371), tolerance = 0.00001)

    # Est
    est <- c(res$res$default$FE[c(4,8,12), "estimate"])
    y <- mean(est)
    expect_equal(tmp$summary$default$FE[4, "M_est"], y, tolerance = 0.001)

    # SE
    se <- c(res$res$default$FE[c(4,8,12), "se"])
    y <- mean(se)
    expect_equal(tmp$summary$default$FE[4, "M_se"], y, tolerance = 0.001)

    # Emperical SE
    expect_equal(tmp$summary$default$FE[4,"SD_est"], sd(est), tolerance = 0.001)


    # Subject_slope
    est <- res$res$default$RE
    est <- est[est$parameter == "subject_slope", "vcov"]

    expect_equal(tmp$summary$default$RE[2, "M_est"], mean(est), tolerance = 0.001)
    expect_equal(tmp$summary$default$RE[2, "prop_zero"], mean(abs(est - 0) < 0.00001),
                 tolerance = 0.001)

    # Subject intercept_slope cor
    est <- res$res$default$RE
    est <- est[est$parameter == "cor_subject", "vcov"]
    expect_equal(tmp$summary$default$RE[5, "M_est"], mean(est), tolerance = 0.001)


    # Cluster slope
    est <- res$res$default$RE
    est <- est[est$parameter == "cluster_slope", "vcov"]

    expect_equal(tmp$summary$default$RE[3, "M_est"], mean(est), tolerance = 0.001)
    expect_equal(tmp$summary$default$RE[3, "prop_zero"], mean(abs(est - 0) < 0.00001),
                 tolerance = 0.001)

    # df_bw
    df <- res$res$default$FE$df_bw
    expect_equal(df[!is.na(df)], c(8,8,8))
})

test_that("simulate with NA para", {
    set.seed(45446)
    p <- study_parameters(n1 = 3,
                          n2 = 5,
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = NA,
                          var_ratio = 0.02,
                          icc_slope = 0.05,
                          cohend = 0
    )


    res <- simulate(p,
                    nsim = 2,
                    formula = sim_formula("y ~ time * treatment + (1 | subject) + (1 | cluster)"))

    res
    x <- summary(res)
    x <- x$summary$default$RE
    # cluster_intercept is NA, but still in moddel formula
    # theta should be 0
    expect_true(x[x$parameter == "cluster_intercept", "theta"] == 0)

})

test_that("simulation summary alpha", {

    set.seed(5446)
    formula <- sim_formula("y ~ treatment * time + (1  | subject)")
    res <- simulate(p, nsim = 10, formula = formula, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- summary(res)[[1]]$default$FE[, "Power"]

    tmp2 <- summary(res, alpha = 0.5)[[1]]$default$FE[, "Power"]

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
    formula <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)")

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    # df_bw
    df <- res$res$default$FE$df_bw
    expect_equal(df[!is.na(df)], c(3,3,3))

    ##
    formula <- sim_formula("y ~ treatment * time + (1 + time | subject) + (1 + treatment:time | cluster)")

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    ##
    formula <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time:treatment | cluster)")

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

    ##
    formula <- sim_formula("y ~ treatment * time + (1 + time | subject) + (1 + time:treatment | cluster)")

    res <- simulate(p, nsim = 3, formula = formula, satterthwaite = FALSE,
                    progress = FALSE, cores = 1, save = FALSE)
    expect_error(summary(res), NA)

})

test_that("simulation random n2", {

    set.seed(5447)
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(func = rnorm(5, 0, 10)),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          cor_subject = -0.5,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = -0.5)
    res <- simulate(p, nsim = 3, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- res$res$default$tot_n[,1]
    expect_gt(length(unique(tmp)), 1)

    # df_bw
    df <- res$res$default$FE$df_bw
    expect_equal(df[!is.na(df)], c(4,4,4))
    expect_error(summary(res), NA)

})


test_that("simulation random n2 some zero", {

    set.seed(5446)
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(func = rnorm(10, 0, 10), replace = 0),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          cor_subject = -0.5,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = -0.5)
    res <- simulate(p, nsim = 3, satterthwaite = FALSE,
                    progress = FALSE)
    tmp <- res$res$default$tot_n[,1]
    expect_length(unique(tmp), 3)

    # df_bw
    df <- res$res$default$FE$df_bw
    expect_gt(length(unique(df[!is.na(df)])), 1)
    expect_error(summary(res), NA)

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
    f <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    res <- simulate(p, nsim = 2, formula = f, satterthwaite = FALSE,
                    progress = FALSE)

    expect_is(res$res$default$FE[1,2], "numeric")
})

# Test simulations run without error with CI
test_that("Simulation runs, unequal_clusters", {
    p <- study_parameters(n1 = 3,
                          n2 = 10,
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)
    f <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)", test = "treatment:time")
    res <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE, CI = TRUE,
                    progress = FALSE)
    expect_is(res$res$default$FE[4, "CI_lwr"], "numeric")
    expect_error(summary(res), NA)
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
    f <- sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    res <- simulate(p, nsim = 2, formula = f, satterthwaite = TRUE,
                    progress = FALSE)

    expect_is(res$res$default$FE[1,2], "numeric")

})






