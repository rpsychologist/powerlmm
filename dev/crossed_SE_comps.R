
p <- study_parameters(design = des,
                      n1 = 3,
                      n2 = 2,
                      n3 = 2,
                      T_end = 2,
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 2, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 1, # time:treatment
                      cor_cluster_intercept_slope = 0,
                      cor_cluster_intercept_intercept_tx = 0,
                      cor_cluster_intercept_slope_tx = 0,
                      cor_cluster_slope_intercept_tx = 0,
                      cor_cluster_slope_slope_tx = 0,
                      cor_cluster_intercept_tx_slope_tx = 0,
                      sigma_error = 10,
                      cor_subject = 0,
                      effect_size = -1
)


d <- simulate_data(p)
fit <- lmer(y ~ time*treatment + (1 + time | subject) + (1 + time + treatment + time:treatment | cluster), data = d)

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

for(i in seq_along(p)) {
   x <- p[[i]]
   p[[i]] <- ifelse(length(x) == 0, 0, x)
}

prepped <- prepare_paras(p)

d <- simulate_data(prepped)
f <- lme4::lFormula(formula = "y ~ time*treatment + (1 + time | subject) + (1 + time + treatment + time:treatment | cluster)",
                    data = d)

pc <- setup_power_calc(object = p,
                       d = d,
                       f = f)
pars <- pc$pars
X <- pc$X
Zt <- pc$Zt
L0 <- pc$L0
Lambdat <- pc$Lambdat
Lind <- pc$Lind

varb <- varb_func(para = pars,
                  X = X,
                  Zt = Zt,
                  L0 = L0,
                  Lambdat = Lambdat,
                  Lind = Lind,
                  crossed = TRUE)
Phi <- varb(Lc = diag(4))

sqrt(Phi)

sqrt(vcov(fit))

# closed
x <- get_pars_short_name(p)
sx <- var_T(p$n1, p$T_end)

varb2 <- with(p, 2*(sigma_error^2 +n1 * sigma_subject_slope^2*sx + 1/2 * n1 * n2 * sigma_cluster_slope_tx^2*sx)/(n1*n2*n3*sx) )

c(sqrt(Phi[4,4]), sqrt(varb2), sqrt(vcov(fit)[4,4]))


