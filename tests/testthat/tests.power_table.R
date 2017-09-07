test_that("power_table", {
    # two level
    # classic se
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              cohend = -0.8)

    # increase only n2
    x <- get_power_table(paras, n2 = 10:15)
    expect_equal(nrow(x), length(10:15))
    expect_equal(colnames(x)[1:2], c("n2", "power"))
    expect_equal(x$tot_n, 10:15*6)
    # expect power to be increasing
    expect_true(all(diff(x$power) > 0))
    expect_is(plot(x), "ggplot")

    # two-level
    # matrix se
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              dropout = dropout_weibull(0.2, 2),
                              cohend = -0.8)

    # increase only n2
    x <- get_power_table(paras, n2 = 10:15)
    expect_equal(nrow(x), length(10:15)*2)
    expect_equal(colnames(x)[1:2], c("n2", "power"))


    x_miss <- x[x$dropout == "with missing", ]
    x_c <- x[x$dropout == "no missing", ]

    expect_equal(x_miss$tot_n, 10:15*6)
    expect_equal(x_c$tot_n, 10:15*6)

    # expect power to be increasing
    expect_true(all(diff(x_miss$power) > 0))
    expect_true(all(diff(x_c$power) > 0))

    # expect complete sample to have higher power
    expect_true(all(x_c$power > x_miss$power))
    expect_is(plot(x), "ggplot")
})

test_that("power_table 2 pars", {
    # two level
    # classic se
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              cohend = -0.8)

    # increase only n2
    x <- get_power_table(paras, n2 = 10:15)
    expect_equal(nrow(x), length(10:15))
    expect_equal(colnames(x)[1:2], c("n2", "power"))

    expect_equal(x$tot_n, 10:15*6)
    # expect power to be increasing
    expect_true(all(diff(x$power) > 0))
    expect_is(plot(x), "ggplot")

    # two-level
    # matrix se
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              dropout = dropout_weibull(0.2, 2),
                              cohend = -0.8)

    # increase n2 and n3
    x <- get_power_table(paras, n2 = 10:11, n3 = 6:8)
    expect_equal(nrow(x), 2*3*2)
    expect_equal(colnames(x)[1:3], c("n2", "n3", "power"))


    x_miss <- x[x$dropout == "with missing", ]
    x_c <- x[x$dropout == "no missing", ]

    # expect complete sample to have higher power
    expect_true(all(x_c$power > x_miss$power))
    expect_is(plot(x), "ggplot")
})


test_that("power_table 3 pars", {
    # two level
    # classic se
    paras <- study_parameters(n1 = 11,
                              n2 = 10,
                              n3 = 6,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              cohend = -0.8)

    # increase only n2
    x <- get_power_table(paras, n2 = 10:11, n3 = 3:4, icc_slope = c(0.1, 0.2))
    expect_equal(nrow(x), 2 * 2 * 2)
    expect_equal(colnames(x)[1:4], c("n2", "n3", "icc_slope", "power"))


    x_miss <- x[x$dropout == "with missing", ]
    x_c <- x[x$dropout == "no missing", ]

    # expect complete sample to have higher power
    expect_true(all(x_c$power > x_miss$power))
    expect_is(plot(x), "ggplot")

    expect_error(get_power_table(paras, n2 = 2:3, n3 = c(2,3), sigma_subject_slope = c(0, 0.2)))
})


# partially nested
test_that("power_table partially nested", {

paras <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 5,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          dropout = dropout_weibull(0.3, 2),
                          cohend = -0.8)
x <- get_power_table(paras, n2 = 2:5, partially_nested = c(TRUE, FALSE))

x_part <- x[x$partially_nested == TRUE, ]
x_full <- x[x$partially_nested == FALSE, ]

# expect partially nested to have more power
expect_true(all(x_part$power > x_full$power))

## standard SE
paras <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = 5,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          cohend = -0.8)
x <- get_power_table(paras, n2 = 2:5, partially_nested = c(TRUE, FALSE))

x_part <- x[x$partially_nested == TRUE, ]
x_full <- x[x$partially_nested == FALSE, ]

# expect partially nested to have more power
expect_true(all(x_part$power > x_full$power))
})

