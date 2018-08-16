library(ggplot2)
library(dplyr)

d <-
    simulate_3lvl_data_crossed(
        n1 = 11,
        n2 = 4,
        n3 = 100,
        T_end = 10,
        fixed_intercept = 4,
        fixed_tx = 0,
        fixed_slope = -1,
        fixed_slope_time_tx = -0.5,
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
        cor_cluster_intercept_tx_slope_tx = 0.5
    )


d <- d %>%
    mutate(y_clust = intercept_cluster + slope_cluster*time + intercept_cluster_tx * treatment + slope_cluster_time_tx * treatment*time)

ggplot(d, aes(time, y_clust, group = subject, color = factor(cluster), linetype = factor(treatment))) + geom_line() + geom_line(aes()) + facet_wrap(~treatment)


ggplot(d, aes(slope_cluster, slope_cluster_time_tx)) + geom_point()

#
ggplot(d, aes(slope_cluster, slope_cluster_time_tx + slope_cluster)) + geom_point()


v1 <- 2
v3 <- 1
cor_cluster_V1V3 <- 0.9
# ICC slopes within clusters
# corr(v1, v1 + v3)
(v1^2 + v1*v3*cor_cluster_V1V3) / (sqrt(v1^2) * sqrt(v1^2 + v3^2 + 2*v1*v3*cor_cluster_V1V3))

#
# corr(v1, v3)


