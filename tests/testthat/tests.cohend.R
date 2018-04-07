## Unfinished file
##
## Draft tests

p0 <- study_parameters(5, 2, 2,
                     T_end = NULL,
                     icc_pre_subject = .5,
                     icc_pre_cluster = 0,
                     effect_size = cohend(0.5, "pretest_SD")
                     )

p0
p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = 20
)
p1
get_power(p1)
get_power(p0)

get_slope_diff(p0)
get_slope_diff(p1)

# how do multiple ES?
p0 <- study_parameters(5, 2, 2,
                 T_end = NULL,
                 icc_pre_subject = .5,
                 icc_pre_cluster = 0,
                 effect_size = c(cohend(0.1, "pretest_SD"),
                                 cohend(0.5, "pretest_SD"))
)

p0


p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = c(10, 20)
)
p1
get_power(p1)


# or vectorize function? this deviates from other funcs
# both should prob work
cohend(c(0.5, 0.1))

# and this
c(0.5, cohend(0.5))



p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = cohend(1)
)
p1

p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = cohend(1:3)
)
p1

# error
p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = c(0, cohend(1:3))
)
p1

p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = c(cohend(0), cohend(1:3))
)
p1
p1 <- study_parameters(5, 2, 2,
                       T_end = NULL,
                       icc_pre_subject = .5,
                       icc_pre_cluster = 0,
                       effect_size = c(cohend(1:2), cohend(3:4))
)
p1


study_parameters(5, unequal_clusters(2,3,4),
                 T_end = NULL,
                 icc_pre_subject = .5,
                 icc_pre_cluster = 0,
                 dropout = dropout_weibull(0.5, 1),
                 effect_size = c(cohend(1:2), cohend(3:4))
)


p1 <- study_parameters(5, 2, 2,
                     T_end = NULL,
                     icc_pre_subject = .5,
                     icc_pre_cluster = 0,
                     effect_size = 0.5
)
p1


study_parameters(5, unequal_clusters(c(5,5,5), c(2,3,4)), 2,
                 T_end = NULL,
                 icc_pre_subject = .5,
                 icc_pre_cluster = 0,
                 effect_size = cohend(0.5, "pretest_SD")
)
