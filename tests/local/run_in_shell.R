library(powerlmm)
library(methods)
library(testthat)
p <- study_parameters(n1 = 11,
                      n2 = 5,
                      n3 = 4,
                      T_end = 10,
                      fixed_intercept = 37,
                      fixed_slope = -0.65,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4726944,
                      sigma_cluster_intercept = 0,
                      sigma_cluster_slope = 0.1084435,
                      sigma_error = 2.8,
                      cor_subject = -0.5,
                      cor_cluster = 0,
                      cohend = -0.8)


res <- simulate(p, nsim = 10, CI = FALSE, satterthwaite = FALSE, cores = 2, progress = FALSE)

test_that("run in shell, no progress", {
    expect_error(res, NA)
    expect_error(summary(res), NA)

})


res <- simulate(p, nsim = 10, CI = FALSE, satterthwaite = FALSE, cores = 2, progress = TRUE)

test_that("run in shell, with progress", {
    expect_error(res, NA)
    expect_error(summary(res), NA)

})
