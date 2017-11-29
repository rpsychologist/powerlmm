
# 3 level -----------------------------------------------------------------


## Multi power
# combine scalar and vector
test_that("power", {
    n2 <- c(3, unequal_clusters(2, 3, 4, 4))
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

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(4, 6))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # n2
  #  x <- unlist(tmp$n2)
  #  names(x) <- NULL
  # expect_equal(x, c("3,3,3", "2,3,4,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})

# combine scalar and per_treatment
test_that("power", {

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

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(4, 5))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2 control
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("3,3,3", "2,3,4"))
    #
    # # n2 treatment
    # x <- unlist(tmp$n2_tx)
    # names(x) <- NULL
    # expect_equal(x, c("3,3,3", "2,2,3,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})
# test df when n3 per_treatment
test_that("power", {
    paras <- study_parameters(n1 = 3,
                              n2 = 5:6,
                              n3 = per_treatment(3, 50),
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              dropout = dropout_weibull(0.3, 2),
                              cohend = -0.8)


    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(51, 51))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("5", "6"))

    # n3 tx
    x <- unlist(tmp$n3_tx)
    names(x) <- NULL
    expect_equal(x, c(50,50))

    # n3 cc
    x <- unlist(tmp$n3_cc)
    names(x) <- NULL
    expect_equal(x, c(3, 3))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})



## Test without linear algebra ---------------------------------------------
test_that("power", {
    n2 <- c(3, per_treatment(5,10))
    paras <- study_parameters(n1 = 3,
                              n2 = n2,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    tdf <- unlist(tmp$df)
    expect_identical(tdf, c(4, 4))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("3", "5"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))

    # compare se with matrix se
    tmp2 <- get_se_3lvl_matrix(as.plcp(paras[1,]))
    expect_equal(tmp$se[[1]], tmp2$se)
})

# combine scalar and per_treatment
test_that("power", {

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
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(4, 5))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2 control
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,3,4"))
    #
    # # n2 treatment
    # x <- unlist(tmp$n2_tx)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,2,3,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))

})
# test df when n3 per_treatment
test_that("power", {
    paras <- study_parameters(n1 = 3,
                              n2 = 5:6,
                              n3 = per_treatment(3, 5),
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.05,
                              cohend = -0.8)


    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(6, 6))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2)
    # names(x) <- NULL
    # expect_equal(x, c("5", "6"))

    # n3 treatment
    x <- unlist(tmp$n3_tx)
    names(x) <- NULL
    expect_equal(x, c(5,5))

    # n3 cc
    x <- unlist(tmp$n3_cc)
    names(x) <- NULL
    expect_equal(x, c(3,3))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))

    # compare se with matrix se
    tmp2 <- get_se_3lvl_matrix(as.plcp(paras[1,]))
    expect_equal(tmp$se[[1]], tmp2$se)
})





# 2 level power ---------------------------------------------------------

## Multi power
# combine scalar and vector
test_that("power", {
    n2 <- c(3, unequal_clusters(2, 3, 4, 4))
    paras <- study_parameters(n1 = 3,
                              n2 = n2,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0,
                              dropout = dropout_weibull(0.3, 2),
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(16, 24))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,3,4,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))

    # compare se with matrix se
    tmp2 <- get_se_3lvl_matrix(as.plcp(paras[1,]))
    expect_equal(tmp$se[[1]], tmp2$se)
})

# combine scalar and per_treatment
test_that("power", {

    n2 <- c(3, per_treatment(unequal_clusters(2, 3, 5),
                             unequal_clusters(2, 2, 3, 4)))

    paras <- study_parameters(n1 = 3,
                              n2 = n2,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0,
                              dropout = dropout_weibull(0.3, 2),
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(16, 19))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2 control
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,3,5"))

    # # n2 treatment
    # x <- unlist(tmp$n2_tx)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,2,3,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})
# test df when n3 per_treatment
test_that("power", {
    paras <- study_parameters(n1 = 3,
                              n2 = 5:6,
                              n3 = per_treatment(3, 10),
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              dropout = dropout_weibull(0.3, 2),
                              cohend = -0.8)


    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(63, 76))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2 control
    # x <- unlist(tmp$n2)
    # names(x) <- NULL
    # expect_equal(x, c("5", "6"))

    # n3 tx
    x <- unlist(tmp$n3_tx)
    names(x) <- NULL
    expect_equal(x, c(10, 10))

    # n3 cc
    x <- unlist(tmp$n3_cc)
    names(x) <- NULL
    expect_equal(x, c(3, 3))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})



# Test without linear algebra ---------------------------------------------
test_that("power", {
    n2 <- c(3, unequal_clusters(2, 3, 4, 4))
    paras <- study_parameters(n1 = 3,
                              n2 = n2,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0,
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(16, 24))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,3,4,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})

# combine scalar and per_treatment
test_that("power", {

    n2 <- c(3, per_treatment(unequal_clusters(2, 3, 5),
                             unequal_clusters(2, 2, 3, 4)))

    paras <- study_parameters(n1 = 3,
                              n2 = n2,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0,
                              cohend = -0.8)

    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(16, 19))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2 control
    # x <- unlist(tmp$n2_cc)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,3,5"))
    #
    # # n2 treatment
    # x <- unlist(tmp$n2_tx)
    # names(x) <- NULL
    # expect_equal(x, c("3", "2,2,3,4"))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})
# test df when n3 per_treatment
test_that("power", {
    paras <- study_parameters(n1 = 3,
                              n2 = 5:6,
                              n3 = per_treatment(3, 10),
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              cohend = -0.8)


    # df
    tmp <- get_power(paras)
    df <- unlist(tmp$df)
    expect_identical(df, c(63, 76))

    # nrow object
    expect_equal(nrow(tmp), 2)

    # # n2
    # x <- unlist(tmp$n2)
    # names(x) <- NULL
    # expect_equal(x, c("5", "6"))

    # n3 treatment
    x <- unlist(tmp$n3_tx)
    names(x) <- NULL
    expect_equal(x, c(10, 10))

    # n3 cc
    x <- unlist(tmp$n3_cc)
    names(x) <- NULL
    expect_equal(x, c(3, 3))

    # set 2 should have more power
    expect_gt(unlist(tmp$power[2]),
              unlist(tmp$power[1]))
})



# partially nested --------------------------------------------------------
test_that("power partially nested", {
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

    p <- update(paras, dropout = 0, partially_nested = TRUE)
    x <- get_power(p)
    x_m <- get_se_3lvl_matrix(p)
    # compare class se and matrix se
    expect_equal(x$se, x_m$se)


    # Test df when per_treatment
    paras <- study_parameters(n1 = 11,
                              n2 = per_treatment(2, 12),
                              n3 = per_treatment(5, 10),
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.1,
                              dropout = dropout_weibull(0.3, 2),
                              cohend = -0.8)
    x <- get_power(paras)
    expect_equal(x$df, 15-2)

})



## New power func
test_varb <- function(object) {
    d <- simulate_data(object)
    f <- lme4::lFormula(formula = create_lmer_formula(object),
                        data = d)

    pc <- setup_power_calc(d, f, object)
    X <- pc$X
    Zt <- pc$Zt
    L0 <- pc$L0
    Lambdat <- pc$Lambdat
    Lind <- pc$Lind

    varb <- varb_func(para = pc$pars, X = X, Zt = Zt, L0 = L0, Lambdat = Lambdat, Lind = Lind)
    bint <- as.numeric(varb(Lc = c(0,0,0,1)))

    old_bint <- get_se_3lvl_matrix(object)$se^2
    expect_equal(old_bint, bint)
}

test_that("varb", {
    p <- study_parameters(n1 = 11,
                              n2 = 4,
                              n3 = 3,
                              T_end = 10,
                              icc_pre_subject = 0.5,
                              icc_pre_cluster = 0,
                              var_ratio = 0.03,
                              icc_slope = 0.1,
                              dropout = 0,
                              cohend = 0.8)

    test_varb(p)

    # unequal
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(2,4,10,50),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          dropout = 0,
                          cohend = 0.8)
    test_varb(p)

    # partially nested
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(2,4,10,50),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE,
                          cohend = 0.8)
    test_varb(p)

    #
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(2,4,10,50),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0.1,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE,
                          cor_cluster = -0.2,
                          cohend = 0.8)
    test_varb(p)


    # Two-level
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(2,4,10,50),
                          n3 = 3,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0,
                          cohend = 0.8)
    test_varb(p)


})



