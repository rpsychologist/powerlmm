p1 <- study_parameters(n1 = n1,
                       n2 = 238/2,
                       T_end = 10,
                       icc_pre_subject = 0.5,
                       cor_subject = -0.5,
                       var_ratio = 0.03,
                       effect_size = cohend(0.5, "posttest_SD"))

p1
p2 <- study_parameters(n1 = n1,
                      n2 = 238/2,
                      T_end = 10,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      var_ratio = 0.03,
                       effect_size = cohend(0.5, "pretest_SD"))

get_power(p1)
get_power(p2)

plot(p1)
plot(p2)

get_slope_diff(p1)
get_slope_diff(p2)


get_sds(p2)



### Slope_sd

D <- 5
n1 <- 10
p <- study_parameters(n1 = n1,
                      n2 = 238/2,
                      T_end = D,
                      sigma_subject_intercept = 0.2,
                      sigma_subject_slope = sqrt(0.0030),
                      sigma_error = sqrt(0.0262),
                      effect_size = sqrt(0.0030) * 0.4 * D)
p
get_power(p)


p1 <- study_parameters(n1 = n1,
                       n2 = 238/2,
                       T_end = D,
                       sigma_subject_intercept = 4,
                       sigma_subject_slope = sqrt(0.0030),
                       sigma_error = sqrt(0.0262),
                       effect_size = cohend(0.4, "slope_SD"))

get_power(p1)


p1 <- study_parameters(n1 = n1,
                       n2 = 238/2,
                       n3 = 4,
                       T_end = D,
                       icc_pre_subject = 0.5,
                       icc_slope = 0.1,
                       var_ratio = 0.05,
                       sigma_error = 10,
                       effect_size = cohend(0.4, "slope_SD"))

p1
p2 <- update(p1, partially_nested = TRUE)

get_slope_diff(p1)
get_slope_diff(p2)

get_power(p1)
get_power(p2)


# subject-levle slope SD
study_parameters(n1 = n1,
                       n2 = 238/2,
                       n3 = 4,
                       T_end = D,
                       icc_pre_subject = 0.5,
                       icc_slope = 0.1,
                       var_ratio = 0.05,
                       sigma_error = 10,
                       effect_size = cohend(0.4, "slope_SD_subject"))
