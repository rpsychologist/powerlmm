
# Errors
test_that("dropout helpers", {
    expect_error(dropout_weibull(1, 1), "'proportion' must be less than 1")
    expect_error(dropout_weibull(-1, 1), "'proportion' can't be negative")
    expect_error(dropout_weibull(0.2, 0), "'rate' must be greater than zero")
    expect_error(dropout_weibull(0.2, -1), "'rate' must be greater than zero")
    expect_is(dropout_weibull(0.2, 2), "plcp_weibull")

    # manual
    expect_error(dropout_manual(0, 0, 1), "Values must be less than 1")
    expect_error(dropout_manual(0, 0, -1), "Values can't be negative")
    expect_error(dropout_manual(0.1, 0.3, 0.5), "Dropout must be 0 at the first time point")
    expect_error(dropout_manual(0, 0.3, 0.2), "Dropout can't be decreasing over time")
    expect_is(dropout_manual(0, 0.2, 0.3)[[1]], "plcp_dropout_manual")

})

test_that("dropout study_parameters", {
    # Non-zero vector
    expect_error(study_parameters(n1 = 10,
                      n2 = 25,
                      icc_pre_subject = 0.5,
                      dropout = 0.5,
                      cohend = 0.5), "'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'")
    expect_error(study_parameters(n1 = 10,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = c(0, 0.5),
                                  cohend = 0.5), "'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'")
    expect_is(study_parameters(n1 = 10,
                     n2 = 25,
                     icc_pre_subject = 0.5,
                     dropout = 0,
                     cohend = 0.5), "plcp")

    # Weibull per_treatment
    expect_error(study_parameters(n1 = 10,
                     n2 = 25,
                     icc_pre_cluster = 0.5,
                     dropout = per_treatment(dropout_weibull(0.2, 2), 0.5),
                     cohend = 0.5), "Treatment group's 'dropout' should be 0 or created by")


    expect_error(study_parameters(n1 = 10,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = per_treatment(0.5, dropout_weibull(0.2, 2)),
                                  cohend = 0.5), "Control group's 'dropout' should be 0 or created by")

    expect_is(study_parameters(n1 = 10,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = per_treatment(0, dropout_weibull(0.2, 2)),
                                  cohend = 0.5), "plcp")

    # Manual per_treatment
    expect_error(study_parameters(n1 = 3,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = per_treatment(0.2, dropout_manual(0, 0.1, 0.3)),
                                  cohend = 0.5), "Control group's 'dropout' should be 0 or created by")

    expect_error(study_parameters(n1 = 3,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = per_treatment(dropout_manual(0, 0.1, 0.3), 0.2),
                                  cohend = 0.5), "Treatment group's 'dropout' should be 0 or created by")


    expect_is(study_parameters(n1 = 3,
                                  n2 = 25,
                                  icc_pre_subject = 0.5,
                                  dropout = per_treatment(dropout_manual(0, 0.1, 0.3), 0),
                                  cohend = 0.5), "plcp")
})


