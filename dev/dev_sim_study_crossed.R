res <- simulate(p,
                formula = sim_formula("y ~ time*treatment + (time | subject) + (time * treatment | cluster)"),
                nsim = 5)

res
summary(res)



d <- simulate_data(p)
fit <- lmer(y ~ time*treatment + (time | subject) + (time * treatment | cluster), data = d)

as.data.frame(VarCorr(fit))



# create lmer syntax
des <- study_design(nested = FALSE)
p <- study_parameters(design = des,
                      n1 = 5,
                      n2 = unequal_clusters(3,4,5),
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
                      cor_cluster_intercept_slope = 0,
                      cor_cluster_intercept_intercept_tx = 0,
                      cor_cluster_intercept_slope_tx = 0,
                      cor_cluster_slope_intercept_tx = 0,
                      cor_cluster_slope_slope_tx = 0,
                      cor_cluster_intercept_tx_slope_tx = 0,
                      effect_size = -1
)


create_lmer_formula(p)



# no cor with cluster_intercept
f <- study_parameters(design = des,
                 n1 = 5,
                 n2 = unequal_clusters(3,4,5),
                 fixed_intercept = 4,
                 fixed_tx = 0,
                 fixed_slope = -1,
                 sigma_subject_intercept = 0,
                 sigma_subject_slope = 0,
                 sigma_cluster_intercept = 1, # cc intercept
                 sigma_cluster_slope = 0, # cc slope
                 sigma_cluster_intercept_tx = 2, # treatment
                 sigma_cluster_slope_tx = 0, # time:treatment
                 sigma_error = 0.1,
                 cor_subject = 0.4,
                 cor_cluster_intercept_slope = NA,
                 cor_cluster_intercept_intercept_tx = NA,
                 cor_cluster_intercept_slope_tx = NA,
                 cor_cluster_slope_intercept_tx = 0,
                 cor_cluster_slope_slope_tx = 0,
                 cor_cluster_intercept_tx_slope_tx = 0,
                 effect_size = -1
) %>%
    create_lmer_formula();f


# no cor with cluster_slope
f <- study_parameters(design = des,
                      n1 = 5,
                      n2 = unequal_clusters(3,4,5),
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 0, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 0, # time:treatment
                      sigma_error = 0.1,
                      cor_subject = 0.4,
                      cor_cluster_intercept_slope = NA,
                      cor_cluster_intercept_intercept_tx = 0,
                      cor_cluster_intercept_slope_tx = 0,
                      cor_cluster_slope_intercept_tx = NA,
                      cor_cluster_slope_slope_tx = NA,
                      cor_cluster_intercept_tx_slope_tx = 0,
                      effect_size = -1
) %>%
    create_lmer_formula();f

lmer(as.formula(f), data = d)

# no cor with cluster_treatment
f <- study_parameters(design = des,
                      n1 = 5,
                      n2 = unequal_clusters(3,4,5),
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 0, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 0, # time:treatment
                      sigma_error = 0.1,
                      cor_subject = 0.4,
                      cor_cluster_intercept_slope = 0,
                      cor_cluster_intercept_intercept_tx = NA,
                      cor_cluster_intercept_slope_tx = 0,
                      cor_cluster_slope_intercept_tx = NA,
                      cor_cluster_slope_slope_tx = 0,
                      cor_cluster_intercept_tx_slope_tx = NA,
                      effect_size = -1
) %>%
    create_lmer_formula();f

lmer(as.formula(f), data = d)

# no cor with cluster_slope_tx
f <- study_parameters(design = des,
                      n1 = 5,
                      n2 = unequal_clusters(3,4,5),
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 0, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 0, # time:treatment
                      sigma_error = 0.1,
                      cor_subject = 0.4,
                      cor_cluster_intercept_slope = 0,
                      cor_cluster_intercept_intercept_tx = 0,
                      cor_cluster_intercept_slope_tx = NA,
                      cor_cluster_slope_intercept_tx = 0,
                      cor_cluster_slope_slope_tx = NA,
                      cor_cluster_intercept_tx_slope_tx = NA,
                      effect_size = -1
) %>%
    create_lmer_formula();f

lmer(as.formula(f), data = d)



f <- study_parameters(design = des,
                      n1 = 5,
                      n2 = unequal_clusters(3,4,5),
                      fixed_intercept = 4,
                      fixed_tx = 0,
                      fixed_slope = -1,
                      sigma_subject_intercept = 0,
                      sigma_subject_slope = 0,
                      sigma_cluster_intercept = 1, # cc intercept
                      sigma_cluster_slope = 0, # cc slope
                      sigma_cluster_intercept_tx = 2, # treatment
                      sigma_cluster_slope_tx = 0, # time:treatment
                      sigma_error = 0.1,
                      cor_subject = 0.4,
                      cor_cluster_intercept_slope = 0,
                      cor_cluster_intercept_intercept_tx = NA,
                      cor_cluster_intercept_slope_tx = NA,
                      cor_cluster_slope_intercept_tx = NA,
                      cor_cluster_slope_slope_tx = NA,
                      cor_cluster_intercept_tx_slope_tx = NA,
                      effect_size = -1
) %>%
    create_lmer_formula();f
