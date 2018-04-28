## single
test_that("stepwise LRT", {
    models <- list(list("label" = "m0",
                        "ll" = -10,
                        "df" = 1),
                   list("label" = "m1",
                        "ll" = -8,
                        "df" = 2),
                   list("label" = "m2",
                        "ll" = -8,
                        "df" = 3),
                   list("label" = "m3",
                        "ll"= -6,
                        "df" = 4)
    )
    expect_equal(step_bw.plcp_sim(models), "m3") # should pick 'm3'
    expect_equal(step_fw.plcp_sim(models), "m1") # should pick 'm1'

    ## multiple
    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1,1)),
                   list("label" = "m1",
                        "ll" = c(-10, -8),
                        "df" = c(2,2)),
                   list("label" = "m2",
                        "ll" = c(-10, -10),
                        "df" = c(3,3)),
                   list("label" = "m3",
                        "ll"= c(-10, -10),
                        "df" = c(4,4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m0", "m1")) # m0, m1
    expect_equal(step_fw.plcp_sim(models), c("m0", "m1")) # m0, m1

    #
    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1,1)),
                   list("label" = "m1",
                        "ll" = c(-10, -8),
                        "df" = c(2,2)),
                   list("label" = "m2",
                        "ll" = c(-6, -10),
                        "df" = c(3,3)),
                   list("label" = "m3",
                        "ll"= c(-10, -10),
                        "df" = c(4,4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m2", "m1")) # m2, m1
    expect_equal(step_fw.plcp_sim(models), c("m0", "m1")) # m0, m1

    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1, 1)),
                   list("label" = "m1",
                        "ll" = c(-8, -8),
                        "df" = c(2, 2)),
                   list("label" = "m2",
                        "ll" = c(-6, -6),
                        "df" = c(3, 3)),
                   list("label" = "m3",
                        "ll"= c(-4, -5),
                        "df" = c(4, 4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m3", "m2")) # m3, m2
    expect_equal(step_fw.plcp_sim(models), c("m3", "m2")) # m3, m2
})


test_that("sim LRT", {
    ## sim
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.5,
                          var_ratio = 0.03)

    f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
    f1 <- sim_formula("y ~ time * treatment + (1 + time || subject)")
    f2 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
    f <- sim_formula_compare("m0" = f0, "m1" = f1, "m2" = f2)


    res <- simulate(p, formula = f, nsim = 4, satterthwaite = FALSE, cores = 1, CI = FALSE)

    x <-  summary(res)
    expect_equal(names(x$summary), c("m0", "m1", "m2"))

    x <- summary(res, model_selection = "FW")
    expect_equal(names(x$summary), "model_selection")
    expect_equal(x$model_direction, "FW")
})

test_that("LRT calcs", {
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.5,
                          var_ratio = 0.05)

    d <- simulate_data(p)
    fit0 <- lme4::lmer(y ~ time * treatment + (1 | subject), data = d)
    fit1 <- lme4::lmer(y ~ time * treatment + (1 + time | subject), data = d)

    av <- anova(fit0, fit1, refit = FALSE)

    ll0 <- stats::logLik(fit0)
    df0 <- attr(ll0, "df")
    ll1 <- stats::logLik(fit1)
    df1 <- attr(ll1, "df")

    m0 <- list("label" = "m0",
               "ll" = as.numeric(ll0),
               "df" = df0)

    m1 <- list("label" = "m1",
               "ll" = as.numeric(ll1),
               "df" = df1)

    expect_equivalent(comp_LRT(m0, m1), av$`Pr(>Chisq)`[2])
})



