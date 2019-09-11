

post_test <- function(fit, d = NULL) {

    # test slopes
    L <- c(rep(1/4, 4), rep(-1/4, 4))
    slope <- lmerTest::contest1D(fit, L = c(rep(0, 8), L))
    slope <- data.frame(parameter = "TE_slope_diff",
                        estimate = slope$Estimate,
                        se = slope$`Std. Error`,
                        pval = slope$`Pr(>|t|)`,
                        df = slope$df,
                        df_bw = NA)

    # CR2
    # fit1 <- tryCatch(nlme::lme(y ~ time * treatment,
    #             random = list(cluster = nlme::pdDiag(~time),
    #                           subject = ~time),
    #             data = d), error = function(e) NA)
    #
    fit1 <- tryCatch(nlme::lme(y ~ time * treatment,
                random = list(subject = ~1),
                data = d), error = function(e) NA)

    cr2 <- tryCatch(clubSandwich::coef_test(fit1, vcov = "CR2"),
                    error = function(e) NA)

    if(is.na(fit1)) {
        cr2 <- data.frame(parameter = "CR2",
                          estimate = NA,
                          se = NA,
                          pval = NA,
                          df = NA,
                          df_bw = NA)
    } else {
        cr2 <- data.frame(parameter = "CR2",
                          estimate = cr2$beta,
                          se = cr2$SE,
                          pval = cr2$p_Satt,
                          df = cr2$df,
                          df_bw = NA)
    }


    rbind(slope, cr2)

}


# Uncorrected fixed effects
f <- sim_formula("y ~ 0 + factor(cluster) + time:factor(cluster) + (1 + time | subject)",
                 test = NULL,
                 post_test = post_test)

f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")

p <- study_parameters(design = study_design(),
                      n1 = 3,
                      n2 = 20,
                      n3 = 4,
                      T_end = 10,
                      icc_pre_subject = 0.5,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      partially_nested = FALSE,
                      effect_size = cohend(c(0.5)))

res <- simulate(p,
                nsim = 100,
                formula = sim_formula_compare("FE" = f,
                                              "LMM" = f2),
                satterthwaite = TRUE,
                cores = 16)

summary(res, para = list("FE" = c("TE_slope_diff", "CR2"),
                         "LMM" = "time:treatment")) %>%
    print(add_cols = c("model", "Parameter"))

