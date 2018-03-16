# test that deparase and n2_num gives correct output




p <- study_parameters(n1 = 3:4,
                      n2 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)
x

p <- study_parameters(n1 = 3:4,
                      n2 = 5,
                      n3 = 2,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)
x




p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(10, 20),
                      n3 = 2,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)


p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(10, 20),
                      n3 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)
x$n3

# pn
p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(10, 20),
                      n3 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      partially_nested = TRUE,
                      cohend = 2
)
p
x <- get_power(p)
x$n3



p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(10, 20),
                      n3 = per_treatment(5,10),
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)
x$n3




p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(unequal_clusters(func = rnorm(10, 10)), unequal_clusters(func = rnorm(5, 15))),
                      n3 = 2,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)

x <- get_power(p)
x$n2_tx_lab
x$n2_cc_lab


x <- get_power(p, R = 3)

p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(unequal_clusters(func = rnorm(10, 10)), unequal_clusters(5,5,5,5,5)),
                      n3 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)
x
x$n3


## rounding


p <- study_parameters(n1 = 3:4,
                      n2 = per_treatment(unequal_clusters(func = rnorm(10, 10)), unequal_clusters(func = rnorm(5, 15))),
                      n3 = 2,
                      T_end = 10,
                      icc_pre_subject = c(0.665, 0.666666666666666),
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
x <- get_power(p)


## two-level
p <- study_parameters(n1 = 4,
                      n2 = per_treatment(unequal_clusters(5,5,5,5,4), unequal_clusters(func = rnorm(10, 10))),
                      n3 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      dropout = 0,
                      cohend = 2
)
p
get_power(p)


p <- study_parameters(n1 = 4,
                      n2 = per_treatment(10, 20),
                      n3 = 5,
                      T_end = 10,
                      icc_pre_subject = 0.6,
                      icc_pre_cluster = 0.05,
                      cor_cluster = -0.5,
                      cor_subject = -0.7,
                      icc_slope = 0.1,
                      var_ratio = 0.02,
                      sigma_error = 10,
                      partially_nested = TRUE,
                      dropout = 0,
                      cohend = 2
)
p

get_power(p)

