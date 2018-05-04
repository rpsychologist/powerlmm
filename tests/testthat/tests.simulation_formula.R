# Test formula
test_that("Simulation formula", {


    f <- check_formula(sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"))
    expect_is(f$default, "plcp_sim_formula")
    expect_equal(f$default$formula, "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_equal(f$default$test, "time:treatment")


    # sim_formula_compare
    f <- check_formula(sim_formula_compare("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"))
    expect_is(f, "plcp_compare_sim_formula")
    expect_equal(f$correct$formula, "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_equal(f$correct$test, "time:treatment")


    # version 1
    f <- sim_formula_compare("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "wrong" = "y ~ treatment * time + (1 + time | subject) + (1 + time | cluster)")
    f <- check_formula(f)
    expect_is(f, "plcp_compare_sim_formula")
    expect_equal(f$correct$formula, "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_equal(f$wrong$formula, "y ~ treatment * time + (1 + time | subject) + (1 + time | cluster)")

    # version 2
    f <- sim_formula_compare("correct" = sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"),
                             "wrong" = sim_formula("y ~ treatment * time + (1 + time | subject) + (1 + time | cluster)"))

    f <- check_formula(f)
    expect_is(f, "plcp_compare_sim_formula")
    expect_equal(f$correct$formula, "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_equal(f$wrong$formula, "y ~ treatment * time + (1 + time | subject) + (1 + time | cluster)")

    # same names
    f <- sim_formula_compare("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Formulas should have unique names.")

    # No names
    expect_error(sim_formula_compare("y ~ treatment * time",
                                     "y ~ treatment * time"), "formula\\(s\\) must have a name")

})

