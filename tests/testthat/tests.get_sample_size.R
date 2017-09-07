

test_that("sampple size helpers", {
    # Balanced
    p <- study_parameters(n1 = 11,
                          n2 = 5,
                          n3 = 7,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n <-5*7
    expect_equivalent(get_tot_n(p), c(n, n, 2*n))
    expect_equivalent(get_n3(p), c(7, 7, 14))
    expect_equivalent(get_n2(p)$treatment, 5)
    expect_equivalent(get_n2(p)$control, 5)

    # unequal clusters
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(15, 20, 33),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n <- 15+20+33
    expect_equivalent(get_tot_n(p), c(n, n, 2*n))
    expect_equivalent(get_n3(p), c(3,3,6))
    expect_equivalent(get_n2(p)$treatment, c(15, 20, 33))
    expect_equivalent(get_n2(p)$control, c(15, 20, 33))

    # per_treatment n2
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(2, 10),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 3*10
    n_cc <- 2*3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(3,3,6))
    expect_equivalent(get_n2(p)$treatment, 10)
    expect_equivalent(get_n2(p)$control, 2)

    # per_treatment n2, unequal_clusters control
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(2, 3), 10),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 3*10
    n_cc <- 2+3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(3, 2, 5))
    expect_equivalent(get_n2(p)$treatment, 10)
    expect_equivalent(get_n2(p)$control, c(2,3))

    # per_treatment n2, unequal_clusters treatment
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(10, unequal_clusters(3, 5)),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 3+5
    n_cc <- 10*3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(2, 3, 5))
    expect_equivalent(get_n2(p)$treatment, c(3,5))
    expect_equivalent(get_n2(p)$control, 10)

    # per_treatment n2, unequal_clusters treatment and control
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(10, 20, 30),
                                             unequal_clusters(3, 5)),
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 3+5
    n_cc <- 10+20+30
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(2, 3, 5))
    expect_equivalent(get_n2(p)$treatment, c(3, 5))
    expect_equivalent(get_n2(p)$control, c(10,20,30))

    # per_treatment n3 only
    p <- study_parameters(n1 = 11,
                          n2 = 5,
                          n3 = per_treatment(2, 10),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 10*5
    n_cc <- 2*5
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(10,2,12))
    expect_equivalent(get_n2(p)$treatment, 5)
    expect_equivalent(get_n2(p)$control, 5)

    # per_treatment n2 and n3
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(7, 6),
                          n3 = per_treatment(2, 10),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1)
    n_tx <- 10*6
    n_cc <- 2*7
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(10,2,12))
    expect_equivalent(get_n2(p)$treatment, 6)
    expect_equivalent(get_n2(p)$control, 7)
})



# partially nested --------------------------------------------------------
test_that("partially nested", {
    # Balanced
    p <- study_parameters(n1 = 11,
                          n2 = 5,
                          n3 = 7,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n <-5*7
    expect_equivalent(get_tot_n(p), c(n, n, 2*n))
    expect_equivalent(get_n3(p), c(7, 0, 7))
    expect_equivalent(get_n2(p)$treatment, 5)
    expect_equivalent(get_n2(p)$control, 5*7)

    # unequal clusters
    p <- study_parameters(n1 = 11,
                          n2 = unequal_clusters(15, 20, 33),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n <- 15+20+33
    expect_equivalent(get_tot_n(p), c(n, n, 2*n))
    expect_equivalent(get_n3(p), c(3,0,3))
    expect_equivalent(get_n2(p)$treatment, c(15, 20, 33))
    expect_equivalent(get_n2(p)$control, c(15 + 20 +33))

    # per_treatment n2
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(2, 10),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 3*10
    n_cc <- 2*3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(3,0,3))
    expect_equivalent(get_n2(p)$treatment, 10)
    expect_equivalent(get_n2(p)$control, 2*3)

    # per_treatment n2, unequal_clusters control
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(2, 3), 10),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 3*10
    n_cc <- 2+3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(3, 0, 3))
    expect_equivalent(get_n2(p)$treatment, 10)
    expect_equivalent(get_n2(p)$control, c(2+3))

    # per_treatment n2, unequal_clusters treatment
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(10, unequal_clusters(3, 5)),
                          n3 = 3,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 3+5
    n_cc <- 10*3
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(2, 0, 2))
    expect_equivalent(get_n2(p)$treatment, c(3,5))
    expect_equivalent(get_n2(p)$control, 10*3)

    # per_treatment n2, unequal_clusters treatment and control
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(unequal_clusters(10, 20, 30),
                                             unequal_clusters(3, 5)),
                          n3 = 5,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 3+5
    n_cc <- 10+20+30
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(2, 0, 2))
    expect_equivalent(get_n2(p)$treatment, c(3, 5))
    expect_equivalent(get_n2(p)$control, c(10+20+30))

    # per_treatment n3 only
    p <- study_parameters(n1 = 11,
                          n2 = 5,
                          n3 = per_treatment(2, 10),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 10*5
    n_cc <- 2*5
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(10, 0, 10))
    expect_equivalent(get_n2(p)$treatment, 5)
    expect_equivalent(get_n2(p)$control, 2*5)

    # per_treatment n2 and n3
    p <- study_parameters(n1 = 11,
                          n2 = per_treatment(7, 6),
                          n3 = per_treatment(2, 10),
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0.03,
                          icc_slope = 0.1,
                          partially_nested = TRUE)
    n_tx <- 10*6
    n_cc <- 2*7
    expect_equivalent(get_tot_n(p), c(n_tx, n_cc, n_tx + n_cc))
    expect_equivalent(get_n3(p), c(10,0,10))
    expect_equivalent(get_n2(p)$treatment, 6)
    expect_equivalent(get_n2(p)$control, 7*2)
})


# Two-level ---------------------------------------------------------------


