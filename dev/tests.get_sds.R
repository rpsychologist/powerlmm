# new tests for partially nested sds


p <- study_parameters(n1 = 11,
                      n2 = per_treatment(20, 40),
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0.25,
                      var_ratio = 0.04,
                      icc_slope = 0.5,
                      partially_nested = TRUE)

get_sds(p, group = "treatment")
get_sds(p, group = "control")


p_tx <- update(p, effect_size = cohend(2, standardizer = "pretest_SD"))
p_cc <- update(p, effect_size = cohend(2, standardizer = "pretest_SD", group = "treatment"))
get_slope_diff(p_tx)
get_slope_diff(p_cc)

get_power(p_tx)
get_power(p_cc)
