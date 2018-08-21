library(lme4)
library(dplyr)

des <- study_design(nested = FALSE)
p <- study_parameters(design = des,
                      n1 = 3,
                      n2 = unequal_clusters(5,5,15,50),
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
                      cor_cluster_slope_intercept_tx = 0.5,
                      cor_cluster_slope_slope_tx = 0.5,
                      cor_cluster_intercept_tx_slope_tx = 0.5,
                      sigma_error = 10,
                      cor_subject = 0,
                      effect_size = -1
)


d <- simulate_data(p)
f <- "y ~ time*treatment + (1 + time | subject) + (1 + time + treatment + treatment:time | cluster)"
fit <- lmer(f, data = d)

vv <- as.data.frame(VarCorr(fit))

get_lme4_sd <- function(x, grp, var1, var2 = NA) {
    x$grp <- gsub("^subject.*", "subject", x$grp)
    x$grp <- gsub("^cluster.*", "cluster", x$grp)

    if(is.na(var2)) {
        sdcor <- dplyr::filter(x, grp == !!grp & var1 == !!var1 & is.na(var2))$sdcor
    } else {
        sdcor <- dplyr::filter(x, grp == !!grp & var1 == !!var1 & var2 == !!var2)$sdcor
    }
    sdcor <- ifelse(length(sdcor) == 0, NA, sdcor)
    ifelse(is.nan(sdcor), 0, sdcor)
}

vv <- as.data.frame(VarCorr(fit))
p$sigma_subject_intercept <- get_lme4_sd(vv, grp = "subject", var1 = "(Intercept)")
p$sigma_subject_slope <- get_lme4_sd(vv, grp = "subject" , var1 = "time")
p$cor_subject <- get_lme4_sd(vv, grp = "subject" , var1 = "(Intercept)", var2 = "time")
p$sigma_cluster_intercept <- get_lme4_sd(vv, grp = "cluster" , var1 = "(Intercept)")
p$sigma_cluster_slope <- get_lme4_sd(vv, grp = "cluster" , var1 = "time")
p$sigma_cluster_intercept_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "treatment")
p$sigma_cluster_slope_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "time:treatment")
p$sigma_error <- vv[vv$grp == "Residual" & is.na(vv$var1) & is.na(vv$var2), "sdcor"]
p$cor_cluster_intercept_slope <- get_lme4_sd(vv, grp = "cluster" , var1 = "(Intercept)", var2 = "time")
p$cor_cluster_intercept_intercept_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "(Intercept)", var2 = "treatment")
p$cor_cluster_intercept_slope_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "(Intercept)", var2 = "time:treatment")
p$cor_cluster_slope_intercept_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "time", var2 = "treatment")
p$cor_cluster_slope_slope_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "time", var2 = "time:treatment")
p$cor_cluster_intercept_tx_slope_tx <- get_lme4_sd(vv, grp = "cluster" , var1 = "treatment", var2 = "time:treatment")



prepped <- prepare_paras(p)

d <- simulate_data(prepped)
f <- lme4::lFormula(formula = f,
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



# closed
x <- get_pars_short_name(p)
sx <- var_T(p$n1, p$T_end)


varb2 <- with(p, 2*(sigma_error^2 +n1 * sigma_subject_slope^2*sx + 1/2 * n1 * n2 * sigma_cluster_slope_tx^2*sx)/(n1*n2*n3*sx) )

c(sqrt(Phi[4,4]),  sqrt(vcov(fit)[4,4]))



# Comp theta
data.frame("make_theta" = round(pc$theta, 4),
           "lme4" = round(getME(fit, "theta"), 4))
