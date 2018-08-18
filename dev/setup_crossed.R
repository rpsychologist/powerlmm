# TODO
# - get_slope_diff: calculate correct SDs for crossed design


library(testthat)
library(dplyr)

des <- study_design(nested = FALSE)

p <- study_parameters(design = des,
    n1 = 5,
    n2 = 4,
    n3 = 10,
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
    cor_cluster_intercept_slope = 0.5,
    cor_cluster_intercept_intercept_tx = 0.5,
    cor_cluster_intercept_slope_tx = 0.5,
    cor_cluster_slope_intercept_tx = 0.5,
    cor_cluster_slope_slope_tx = 0.9,
    cor_cluster_intercept_tx_slope_tx = 0.5,
    effect_size = -1
)

get_slope_diff(p)
d <- simulate_data(p)

# tests
expect_length(unique(d$cluster), p$n3)
expect_length(unique(d$subject), p$n3 * p$n2)

## Unequal clusters
### need to have same n3 in both arms

### expect error
p <- study_parameters(design = des,
                      n1 = 5,
                      n2 = 4,
                      n3 = per_treatment(4,5),
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
                      effect_size = -1
)


## Unequal n2 per treatment
### i.e. therapists with unbalanced patients in tx A or B

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
                      effect_size = -1
)
d <- simulate_data(p)


d %>% group_by(treatment,cluster) %>%
    summarise(n = length(unique(subject)))


p <- study_parameters(design = des,
                      n1 = 5,
                      n2 = per_treatment(control = unequal_clusters(3, 4, 5),
                                         treatment = unequal_clusters(10, 10, 10)),
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
                      effect_size = -1
)
d <- simulate_data(p)

d %>% group_by(treatment,cluster) %>%
    summarise(n = length(unique(subject)))


## what if not all clusters are crossed?
## different number of clusters per treatment arm
## give error?
p <- study_parameters(design = des,
                      n1 = 5,
                      n2 = per_treatment(unequal_clusters(3,4,5),
                                         unequal_clusters(10,10,10, 10)),
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
                      effect_size = -1
)

simulate_data(p)

## Missing data
p <- study_parameters(design = des,
                      n1 = 11,
                      n2 = per_treatment(control = unequal_clusters(30, 40, 50),
                                         treatment = unequal_clusters(90, 110, 120)),
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
                      dropout = dropout_weibull(0.4, 1),
                      effect_size = -1
)

d <- simulate_data(p)

d %>%
    filter(time == 10) %>%
    group_by(treatment) %>%
    summarise(mean(is.na(y)))

d %>%
    filter(time == 10) %>%
    group_by(cluster, treatment) %>%
    summarise(mean(is.na(y)))
