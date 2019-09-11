library(brms)

nsim <- 5
cores <- 2
# sim parameters ----------------------------------------------------------
des <- study_design(family = "hurdle")

p <- study_parameters(
    design = des,
    n1 = 5,
    n2 = 50,
    T_end = 10,
    fixed_intercept = log(30), # median(Y > 0)
    fixed_hu_intercept = qlogis(0.8), # prop == 0
    fixed_slope = log(0.99),
    fixed_hu_slope = log(1),
    sd_hu_intercept = 3,
    sd_hu_slope = 0.2,
    sd_intercept = 2,
    sd_slope = 0.05,
    cor_intercept_slope = -0.15,
    cor_intercept_hu_intercept = -0.66,
    cor_intercept_hu_slope = 0.2,
    cor_slope_hu_intercept = -0.1,
    cor_slope_hu_slope = -0.1,
    cor_hu_intercept_hu_slope = 0.15,
    shape = 1.6,
    RR_cont = 0.33,
    OR_hu = 2,
    marginal = TRUE,
    dropout = dropout_weibull(0.4, 1),
    family = "gamma")


# Models ------------------------------------------------------------------
d <- simulate_data(p)
bfit_gamma_ctp <- brm(bf(y ~ time * treatment + (1 + time | c | subject),
                         hu ~ time * treatment + (1 + time | c | subject)),
                      data = d,
                      family = hurdle_gamma,
                      chains = 1,
                      cores = 1,
                      silent = TRUE,
                      iter = 1)


bfit_gamma_mtp <- brm(bf(y ~ time * treatment + (1 + time | c | subject),
                         hu ~ time * treatment + (1 + time | c | subject)),
                      data = d,
                      family = hurdle_gamma_mtp,
                      stanvar = hurdle_gamma_mtp_stanvars,
                      chains = 1,
                      cores = 1,
                      iter = 1)

f0 <- sim_formula(bfit_gamma_mtp, iter = 2000, marginalize = TRUE)


res_gamma <- simulate(p,
                      formula = sim_formula_compare("hurdle_mtp" = f0),
                      nsim = 200,
                      cores = 16)

summary(res_gamma, marginalize = TRUE)

saveRDS(res_gamma, "sim_gamma_hurdle.rds")

