# TODO
# - get_slope_diff: calculate correct SDs for crossed design



des <- study_design(nested = FALSE)

p <- study_parameters(design = des,
    n1 = 11,
    n2 = 4,
    n3 = 100,
    T_end = 10,
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

## Unequal clusters
### need to have same n3 in both arms

## Unequal n2 per treatment
### i.e. therapists with unbalanced patients in tx A or B

## Missing data


