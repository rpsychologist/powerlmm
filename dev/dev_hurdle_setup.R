des <- structure(list(), class = "plcp_hurdle")

p <- study_parameters(
    design = des,
    n1 = 10,
    n2 = per_treatment(10, 50),
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
    dropout = per_treatment(control = dropout_weibull(0.3, 1),
                            treatment = dropout_weibull(0.5, 1)
                            ),
    family = "gamma")

class(p)
create_R_cov(p)

d <- simulate_data(p)

# missing
d %>%
    group_by(treatment, time) %>%
    summarise(mean(is.na(y))) %>%
    print(n=222)

# n2
d %>%
    group_by(treatment, time) %>%
    summarise(n())

get_n2(p)



## Multi
pm <- study_parameters(
    design = des,
    n1 = 10,
    n2 = c(5, 10),
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
    dropout = per_treatment(control = dropout_weibull(0.3, 1),
                            treatment = dropout_weibull(0.5, 1)
    ),
    family = "gamma")

pm
