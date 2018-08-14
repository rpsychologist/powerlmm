# TODO
# * custom print methods

study_parameters()
des_nested <- study_design()
des_crossed <- study_design(nested = FALSE)

study_parameters(des_nested, )
study_parameters(des_crossed, )

p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 6,
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.05,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3,
                                                rate = 1/2),
                      partially_nested = TRUE,
                      effect_size = cohend(-0.8,
                                           standardizer = "pretest_SD"))
