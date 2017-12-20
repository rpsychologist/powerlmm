test_that("get_DEFT", {
    p <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0L,
                          icc_slope = 0.1,
                          sigma_error = 1L,
                          var_ratio = 0.05)

    DEFT1 <- get_DEFT(p)

    se1 <- get_power(p)$se
    se2 <- get_power(update(p, icc_slope = 0))$se

    DEFT2 <- se1/se2

    expect_equal(DEFT1$DEFT, DEFT2)
    expect_error(DEFT1, NA)
})
