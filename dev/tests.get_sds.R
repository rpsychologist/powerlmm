# new tests for partially nested sds


p <- study_parameters(n1 = n1,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0.25,
                      var_ratio = 0.04,
                      icc_slope = 0.5,
                      partially_nested = TRUE)

get_sds(p, group = "treatment")
get_sds(p, group = "control")
