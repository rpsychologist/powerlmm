
p <- study_parameters(design = des,
                      n1 = 3,
                      n2 = 5,
                      n3 = 5,
                      T_end = 2,
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 20,
                      sigma_subject_slope = 2,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 2, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 1, # time:treatment
                      cor_cluster_intercept_slope = 0.5,
                      cor_cluster_intercept_intercept_tx = 0.5,
                      cor_cluster_intercept_slope_tx = 0.5,
                      cor_cluster_slope_intercept_tx = 0.5,
                      cor_cluster_slope_slope_tx = 0.5,
                      cor_cluster_intercept_tx_slope_tx = 0.5,
                      sigma_error = 20,
                      cor_subject = NA,
                      effect_size = -1
)


d <- simulate_data(p)
fit <- lmer(y ~ time*treatment + (1 + time || subject) + (1 + time + treatment + time:treatment | cluster), data = d)

vv <- as.data.frame(VarCorr(fit))

vv$grp <- gsub("^subject.*", "subject", vv$grp)
vv$grp <- gsub("^cluster.*", "cluster", vv$grp)


p$sigma_subject_intercept <- filter(vv, grp == "subject" & var1 == "(Intercept)" & is.na(var2))$sdcor
p$sigma_subject_slope <- filter(vv, grp == "subject" & var1 == "time" & is.na(var2))$sdcor
p$sigma_cluster_intercept <- filter(vv, grp == "cluster" & var1 == "(Intercept)" & is.na(var2))$sdcor
p$sigma_cluster_slope <- filter(vv, grp == "cluster" & var1 == "time" & is.na(var2))$sdcor
p$sigma_cluster_intercept_tx <- filter(vv, grp == "cluster" & var1 == "treatment" & is.na(var2))$sdcor
p$sigma_cluster_slope_tx <- filter(vv, grp == "cluster" & var1 == "time:treatment" & is.na(var2))$sdcor
p$sigma_error <- filter(vv, grp == "Residual" & is.na(var1) & is.na(var2))$sdcor
p$cor_cluster_intercept_slope <- filter(vv, grp == "cluster" & var1 == "(Intercept)" & var2 == "time")$sdcor
p$cor_cluster_intercept_intercept_tx <- filter(vv, grp == "cluster" & var1 == "(Intercept)" & var2 == "treatment")$sdcor
p$cor_cluster_intercept_slope_tx <- filter(vv, grp == "cluster" & var1 == "(Intercept)" & var2 == "time:treatment")$sdcor
p$cor_cluster_slope_intercept_tx <- filter(vv, grp == "cluster" & var1 == "time" & var2 == "treatment")$sdcor
p$cor_cluster_slope_slope_tx <- filter(vv, grp == "cluster" & var1 == "time" & var2 == "time:treatment")$sdcor
p$cor_cluster_intercept_tx_slope_tx <- filter(vv, grp == "cluster" & var1 == "treatment" & var2 == "time:treatment")$sdcor
p$cor_subject <- 0




###
pars <- get_pars_short_name(p)
pars <- with(pars,
             list(u0 = u0^2,
                  u1 = u1^2,
                  u01 = u0 * u1 * u01,
                  v0 = v0^2,
                  v1 = v1^2,
                  v2 = v2^2,
                  v3 = v3^2,
                  v01 = v0 * v1 * v01,
                  v02 = v0 * v2 * v02,
                  v03 = v0 * v3 * v03,
                  v12 = v1 * v2 * v12,
                  v13 = v1 * v3 * v13,
                  v23 = v2 * v3 * v23,
                  sigma = p$sigma_error^2)
)


sigma <- sqrt(pars$sigma)
lvl2 <- make_theta_vec(pars$u0, pars$u01, pars$u1)/sigma

x <- with(pars, matrix(c(v0, v01, v02,  v03,
                      v01, v1,  v12, v13,
                      v02, v12, v2,   v23,
                      v03, v13, v23,  v3), ncol = 4))
x[is.na(x)] <- 0
L <- suppressWarnings(chol(x))


m <- L/sigma

crossprod(L2)
crossprod(m2 * sigma)

round(x, 3)

# Same if positive definite
m <- t(L/sigma)
m[lower.tri(m, diag = TRUE)]
round(m[lower.tri(m, diag= TRUE )], 3)
round(getME(fit, "theta"), 3)[-c(1:2)]

# Pivot
L <- suppressWarnings(chol(x, pivot=TRUE))

pivot <- attr(L, "pivot")
L2 <- L[, order(pivot)]

m2 <- L2/sigma

round(, 3)
round(getME(fit, "theta"), 3)[-c(1:2)]
