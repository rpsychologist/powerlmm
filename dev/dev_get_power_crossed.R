p <- study_parameters(design = study_design(nested = FALSE),
                      n1 = 5,
                      n2 = 5,
                      n3 = 8,
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 2, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 1, # time:treatment
                      sigma_error = 0.1,
                      cor_subject = 0.4,
                      dropout = dropout_weibull(0.3, 1),
                      effect_size = -1
)

get_power(p)

res <-
