p <- study_parameters(design = study_design(nested = FALSE),
                      n1 = 11,
                      n2 = 5,
                      n3 = 12,
                      T_end = 10,
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 10,
                      sigma_subject_slope = 1.4,
                      sigma_cluster_intercept = NA, # cc intercept
                      sigma_cluster_slope = 0.5, # cc slope
                      sigma_cluster_intercept_tx = NA, # treatment
                      sigma_cluster_slope_tx = 0.4, # time:treatment
                      sigma_error = 10,
                      cor_subject = 0.4,
                      effect_size = -5
)

get_power(p)


res <- simulate(p,
                formula = sim_formula("y ~ time*treatment + (1 + time | subject) + (0  + time + time:treatment || cluster)"),
                nsim = 3000, cores = 16, satterthwaite = FALSE)

summary(res, para = "time:treatment")


# SE
x <- get_pars_short_name(p)
sx <- var_T(11, 10)

varb <- with(p, 2*(sigma_error^2 +n1 * sigma_subject_slope^2*sx + 1/2 * n1 * n2 * sigma_cluster_slope_tx^2*sx)/(n1*n2*n3*sx) )
sqrt(varb)


## lmer fit
fit <- lmer(y ~ time*treatment + (1 + time | subject) + (1  + time * treatment | cluster), data = d )

getME(fit, "theta")


x <- as.data.frame(VarCorr(fit))

m <- matrix(c(x$vcov[4], x$vcov[8], x$vcov[9], x$vcov[10],
         x$vcov[8], x$vcov[5], x$vcov[11], x$vcov[12],
         x$vcov[9], x$vcov[11], x$vcov[6], x$vcov[13],
         x$vcov[10], x$vcov[12], x$vcov[13], x$vcov[7]), ncol = 4)

L <- suppressWarnings(chol(m, pivot = TRUE))
L <- L[, order(attr(L, "pivot"))]
L/getME(fit, "sigma")

getME(fit, "theta")[-c(1:3)]
