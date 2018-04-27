compare_sim_formulas("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
                     "wrong" = "y ~ treatment * time + (1 + time | subject)") %>%
    class


compare_sim_formulas("correct" = sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"),
                     "wrong" =  sim_formula("y ~ treatment * time + (1 + time | subject)"))


class(sim_formula("a"))


p <- study_parameters(n1 = 3,
                      n2 = 3,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      var_ratio = 0)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)


# should give error
res <- simulate(p, formula = list("wrong" = "y ~ time * treatment + (1 | subject)",
                                  "corredt" = "y ~ time * treatment + (1 | subject)"), nsim = 1)


res
res <- simulate(p, formula = compare_sim_formulas("correct" = sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"),
                                                  "correct" =  sim_formula("y ~ treatment * time + (1 + time | subject)")),
                nsim = 1)



res <- simulate(p, formula = compare_sim_formulas(sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")),
                nsim = 3)

summary(res)

res <- simulate(p, formula = sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"),
                nsim = 1)

# should work
compare_sim_formulas("a" = sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"))

# should work
compare_sim_formulas("a" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")

# expect error
compare_sim_formulas(sim_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"))
compare_sim_formulas("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")






# old tests ---------------------------------------------------------------
# need updating

# Test formula
test_that("Simulation formula (character)", {

    expect_is(check_formula("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"), "list")
    expect_error(check_formula("y ~ treatment2 * time + (1 + time | subject) + (0 + time | cluster)"))
    expect_is(check_formula("y ~ treatment + time + treatment:time + (1 + time | subject) + (0 + time | cluster)"), "list")
    expect_error(check_formula("y ~ treatment + time + treatment:time + (1 + time2 | subject) + (0 + time | cluster)"))
    expect_error(check_formula("log(y) ~ treatment + time + treatment:time + (1 + time | subject) + (0 + time | cluster)"))
    expect_error(check_formula("log(y) ~ treatment + poly(time, 2) + treatment:time + (1 + time | subject) + (0 + time | cluster)"))

})
# Test formula
test_that("Simulation formula (list)", {


    expect_is(check_formula(list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")), "list")
    expect_is(check_formula(list("wrong"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")), "list")

    f <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "wrong" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_is(check_formula(f), "list")

    # same names
    f <- list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Both formulas can't have the same name")

    # Wrong names
    f <- list("correct"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "Wong"="y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Formula names must be either 'correct' or 'wrong'")

    # No names
    f <- list("y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)",
              "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "Formula should be a named list")

    # Wrong terms
    f <- list("correct" = "y ~ treatment2 * time + (1 + time | subject) + (0 + time | cluster)",
              "wrong" = "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)")
    expect_error(check_formula(f), "treatment2 is not an allowed variable name.")

})

