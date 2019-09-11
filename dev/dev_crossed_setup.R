des <- study_design(nested = TRUE)
p1 <- study_parameters(design = des,
                      n1 = 3:4,
                      n2 = unequal_clusters(5,5,15,50),
                      T_end = 2,
                      fixed_intercept = 4,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_error = 10,
                      cor_subject = 0,
                      effect_size = -1
)

as.plcp(p1[1,])
p1



des <- study_design(nested = FALSE)
p1 <- study_parameters(design = des,
                       n1 = 3:4,
                       n2 = unequal_clusters(5,5,15,50),
                       T_end = 2,
                       fixed_intercept = 4,
                       fixed_slope = -1,
                       sigma_subject_intercept = 0,
                       sigma_error = 10,
                       cor_subject = 0,
                       effect_size = -1
)

as.plcp(p1[1,])
p1



# dropout -----------------------------------------------------------------

des <- study_design(nested = TRUE)
p1 <- study_parameters(design = des,
                       n1 = 10,
                       n2 = unequal_clusters(5,5,15,50),
                       fixed_intercept = 4,
                       fixed_slope = -1,
                       sigma_subject_intercept = 0,
                       sigma_error = 10,
                       cor_subject = 0,
                       dropout = dropout_weibull(0.3, 1),
                       effect_size = -1
)

p1


