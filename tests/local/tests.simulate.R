library(dplyr)
message("Running simulations...")
source("run_simulations.R")
res <- readRDS("simres.rds")
message("Running tests...")
test_local_sim <- function(res, R = 1, cores = 1) {
    x <- summary(res)
    p <- res$p
    pow <- get_power(p, df = "satterth", R = cores, cores = cores)

    RB_RE <- x$summary$correct$RE %>%
        mutate(RB = (M_est - theta)/theta) %>%
        .$RB

    expect_true(all(RB_RE < 0.05))

    # Fixed SE
    RB_fixed_SE <- x$summary$correct$FE %>%
        mutate(RB = (M_se - SD_est)/SD_est) %>%
        .$RB
    expect_true(all(RB_fixed_SE < 0.05))

    # Fixed theta
    slope_diff <- with(p, cohend/T_end * sqrt(sigma_subject_intercept^2 + sigma_error^2))
    theta <- with(p, c(fixed_intercept, fixed_slope, slope_diff))

    RB_fixed_theta <- x$summary$correct$FE[c(1,3,4), ] %>%
        mutate(RB = (M_est - theta)/theta) %>%
        .$RB

    expect_true(all(RB_fixed_theta < 0.05))

    # Compare power
    diff_power <- abs(pow$power - x$summary$correct$FE[4, "Power_satt"])
    expect_lt(diff_power, 0.01)
}

test_that("res1 few_clusters", {
    test_local_sim(res[[1]])
})
test_that("res2 few_clusters PN", {
    test_local_sim(res[[2]])
})
test_that("res3 more_clusters", {
    test_local_sim(res[[3]])
})
test_that("res4 unequal_clusters", {
    test_local_sim(res[[4]])
})
test_that("res5 random_clusters", {
    test_local_sim(res[[5]], R = 100, cores = 30)
})
