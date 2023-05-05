suppressMessages(library(dplyr))
test_that("number of subjects, uneqal", {
    p <- study_parameters(n1 = 5,
                          n2 = unequal_clusters(5, 7, 15),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n
    x <- rep(c(5,7,15), each = 5)
    rep(x, 2)
    expect_equal(n, rep(x,2))
})


test_that("number of subjects, uneqal, dropout", {
    p <- study_parameters(n1 = 5,
                          n2 = unequal_clusters(5, 7, 15),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)


    d <- simulate_data(p)

    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n

    x <- rep(c(5,7,15), each = 5)
    rep(x, 2)

    expect_equal(n, rep(x,2))
})
test_that("number of subjects, per_treatment, uneqal", {
    n2 <- per_treatment(unequal_clusters(5, 7, 15),
                  unequal_clusters(2, 4, 10, 20))
    p <- study_parameters(n1 = 5,
                          n2 = n2,
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)


    d <- simulate_data(p)

    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n

    x <- rep(c(5,7,15), each = 5)
    x2 <-  rep(c(2, 4, 10, 20), each = 5)

    expect_equal(n, c(x, x2))
})
test_that("number of subjects, per_treatment, uneqal, dropout", {
    n2 <- per_treatment(unequal_clusters(5, 7, 15),
                  unequal_clusters(2, 4, 10, 20))
    p <- study_parameters(n1 = 5,
                          n2 = n2,
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n

    x <- rep(c(5,7,15), each = 5)
    x2 <-  rep(c(2, 4, 10, 20), each = 5)

    expect_equal(n, c(x, x2))
})

test_that("number of clusters, per_treatment", {

    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = per_treatment(10, 5),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n

    expect_equal(n, rep(5, 5*10 + 5*5))

    n3 <- d %>%
        group_by(treatment, time) %>%
        summarise(n = length(unique(cluster))) %>%
        .$n
    expect_equal(n3, rep(c(10, 5), each = 5))
})
test_that("number of clusters, per_treatment, dropout", {

    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          n3 = per_treatment(10, 5),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n

    expect_equal(n, rep(5, 5*10 + 5*5))

    n3 <- d %>%
        group_by(treatment, time) %>%
        summarise(n = length(unique(cluster))) %>%
        .$n
    expect_equal(n3, rep(c(10, 5), each = 5))
})

# both n2 and n3 per treatment
test_that("number of clusters, per_treatment", {

    p <- study_parameters(n1 = 5,
                          n2 = per_treatment(2, 5),
                          n3 = per_treatment(10, 5),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, cluster, time) %>%
        summarise(n = n()) %>%
        .$n
    x <- rep(2, 5*10)
    x2 <- rep(5, 5*5)
    expect_equal(n, c(x, x2))

    n3 <- d %>%
        group_by(treatment, time) %>%
        summarise(n = length(unique(cluster))) %>%
        .$n
    expect_equal(n3, rep(c(10, 5), each = 5))
})

# Weibull
test_that("proportion of dropout", {

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout_weibull(0.2, 1),
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[10], 0.2)
    expect_equal(tx$miss[10], 0.2)
})

test_that("proportion of dropout, per_treatment", {

    dropout <- per_treatment(dropout_weibull(0.2, 1),
                            dropout_weibull(0.3, 3))

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[10], 0.2, tolerance = 0.01)
    expect_equal(tx$miss[10], 0.3, tolerance = 0.01)
})
test_that("proportion of dropout, only in control", {

    dropout <- per_treatment(dropout_weibull(0.2, 1),
                             0)

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[10], 0.2, tolerance = 0.01)
    expect_equal(tx$miss[10], 0, tolerance = 0.01)
})
test_that("proportion of dropout, only in control", {

    dropout <- per_treatment(dropout_weibull(0.2, 1),
                             0)

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          deterministic_dropout = FALSE,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)

    expect_equal(tx$miss[10], 0, tolerance = 0.01)
})
test_that("proportion of dropout, only in tx", {

    dropout <- per_treatment(0,
                             dropout_weibull(0.2, 1))

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[10], 0, tolerance = 0.01)
    expect_equal(tx$miss[10], 0.2, tolerance = 0.01)
})
test_that("proportion of dropout, only in tx", {

    dropout <- per_treatment(0,
                             dropout_weibull(0.2, 1))

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          deterministic_dropout = FALSE,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    cc <- n %>% filter(treatment == 0)

    expect_equal(cc$miss[10], 0, tolerance = 0.001)
})

## Manual dropout
test_that("proportion of dropout", {

    dropout <- dropout_manual(seq(0,0.5, length.out = 10))

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(dp$control, seq(0,0.5, length.out = 10))
    expect_equal(dp$treatment, seq(0,0.5, length.out = 10))
    expect_equal(cc$miss[10], 0.5, tolerance = 0.01)
    expect_equal(tx$miss[10], 0.5, tolerance = 0.01)
})

# Manual dropout, per treatment
test_that("proportion of dropout, per_treatment", {

    dropout <- per_treatment(dropout_manual(seq(0,0.5, length.out = 10)),
                             dropout_manual(seq(0,0.3, length.out = 10)))

    p <- study_parameters(n1 = 10,
                          n2 = 25,
                          n3 = 25,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(dp$control, seq(0,0.5, length.out = 10))
    expect_equal(dp$treatment, seq(0,0.3, length.out = 10))
    expect_equal(cc$miss[10], 0.5, tolerance = 0.01)
    expect_equal(tx$miss[10], 0.3, tolerance = 0.01)
})


test_that("proportion of dropout, per_treatment", {

    p <- study_parameters(n1 = 10,
                          n2 = 5,
                          n3 = 5,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = 0.5)

    # partially nested
    set.seed(3)
    d <- simulate_data(p)
    n <- d %>%
        group_by(treatment) %>%
        summarise(n = length(unique(cluster_slope))) %>%
        .$n

    expect_equal(n, c(1,5))

    # not partially_nested
    p <- update(p, partially_nested = FALSE)

    d <- simulate_data(p)
    n <- d %>%
        group_by(treatment) %>%
        summarise(n = length(unique(cluster_slope))) %>%
        .$n
    expect_equal(n, c(5,5))

})



# Multi para
test_that("multi_para #1", {
    p <- study_parameters(n1 = 5:6,
                          n2 = unequal_clusters(5, 7, 15),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          cohend = 0.5)


    d <- simulate_data(p)
    n <- d %>% group_by(subject) %>%
        summarise(n = n()) %>%
        .$n

    expect_true(all(n == 5))
})

# Different dropout patterns
test_that("multi_para #1", {

    dropout <- c(dropout_weibull(0.1, 2),
                 dropout_weibull(0.2, 0.5))
    p <- study_parameters(n1 = 5:6,
                          n2 = unequal_clusters(10, 25, 40, 50, 80, 100, 100),
                          n3 = 2,
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          dropout = dropout,
                          cohend = 0.5)

    set.seed(5444)
    d <- simulate_data(p)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p)

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[5], 0.1, tolerance = 0.01)
    expect_equal(tx$miss[5], 0.1, tolerance = 0.01)

    ## data set #2
    d <- simulate_data(p, n = 2)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    n1 <- d %>%
        group_by(subject) %>%
        summarise(n = n()) %>%
        .$n

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p, n = 2)

    expect_true(all(n1 == 6))

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[6], 0.1, tolerance = 0.01)
    expect_equal(tx$miss[6], 0.1, tolerance = 0.01)

    ## data set #3 diff dropout
    d <- simulate_data(p, 3)
    n <- d %>% group_by(treatment, time) %>%
        summarise(miss = mean(is.na(y)))

    n1 <- d %>%
        group_by(subject) %>%
        summarise(n = n()) %>%
        .$n

    tx <- n %>% filter(treatment == 1)
    cc <- n %>% filter(treatment == 0)

    dp <- get_dropout(p, 3)

    expect_true(all(n1 == 5))

    expect_equal(tx$miss, dp$treatment, tolerance = 0.01)
    expect_equal(cc$miss, dp$control, tolerance = 0.01)
    expect_equal(cc$miss[5], 0.2, tolerance = 0.01)
    expect_equal(tx$miss[5], 0.2, tolerance = 0.01)



})

test_that("proportion of dropout, per_treatment", {

    p <- study_parameters(n1 = 10,
                          n2 = unequal_clusters(func = rpois(5, lambda = 5)),
                          sigma_subject_intercept = 1.44,
                          icc_pre_cluster = 0,
                          sigma_subject_slope = 0.2,
                          icc_slope = 0.05,
                          sigma_error = 1.44,
                          partially_nested = TRUE,
                          cohend = 0.5)

    # partially nested
    set.seed(3)
    d <- simulate_data(p)
    n <- d %>%
        group_by(treatment) %>%
        summarise(n = length(unique(cluster_slope))) %>%
        .$n

    expect_equal(n, c(1,5))

    n2 <- d %>%
        group_by(treatment, cluster) %>%
        filter(time == 0) %>%
        summarise(n = length(cluster)) %>%
        .$n

    expect_true(length(unique(n2)) > 1)

    # not partially_nested
    p <- update(p, partially_nested = FALSE)

    d <- simulate_data(p)
    n <- d %>%
        group_by(treatment) %>%
        summarise(n = length(unique(cluster_slope))) %>%
        .$n
    expect_equal(n, c(5,5))

    n2 <- d %>%
        group_by(treatment, cluster) %>%
        filter(time == 0) %>%
        summarise(n = length(cluster)) %>%
        .$n

    expect_true(length(unique(n2)) > 1)


    # test prepped
    prepped <- prepare_paras(p)

    d <- simulate_data(prepped)
    n <- d %>%
        group_by(treatment) %>%
        summarise(n = length(unique(cluster_slope))) %>%
        .$n
    expect_equal(n, c(5,5))

    n2 <- d %>%
        group_by(treatment, cluster) %>%
        filter(time == 0) %>%
        summarise(n = length(cluster)) %>%
        .$n
    expect_true(length(unique(n2)) > 1)
    n2_2 <- as.integer(get_n2(prepped$control)$treatment)
    expect_identical(rep(n2_2, 2), n2)


})

