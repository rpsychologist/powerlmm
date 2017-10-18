test_that("create_lmer_formula", {
    p <- study_parameters(n1 = 5,
                          n2 = 2,
                          n3 = 2,
                          icc_pre_subject = 0.5,
                          icc_pre_cluster = 0,
                          var_ratio = 0,
                          icc_slope = 0,
                          cohend = 0.5)

    # Two-level
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 | subject)")

    p <- update(p, var_ratio = 0.1)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time || subject)")

    p <- update(p, cor_subject = -0.5)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject)")

    p <- update(p, cor_subject = 0, icc_pre_subject = 0)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (0 + time | subject)")

    # Three-level
    p <- update(p, icc_pre_subject = 0.3, var_ratio = 0, icc_pre_cluster = 0.1)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 | subject) + (1 | cluster)")

    p <- update(p, icc_slope = 0.05, var_ratio = 0.1)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time || subject) + (1 + time || cluster)")

    p <- update(p, cor_cluster = 0.2)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time || subject) + (1 + time | cluster)")


    p <- update(p, cor_subject = 0.3)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject) + (1 + time | cluster)")

    p <- update(p, icc_pre_cluster = 0)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)")

    # Partially nested
    p <- update(p, partially_nested = TRUE)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject) + (0 + treatment:time | cluster)")

    p <- update(p, icc_pre_cluster =  0.1)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject) + (0 + treatment + treatment:time | cluster)")

    p <- update(p, cor_cluster = 0)
    expect_equivalent(create_lmer_formula(p),
                      "y ~ time*treatment + (1 + time | subject) + (0 + treatment + treatment:time || cluster)")

})
