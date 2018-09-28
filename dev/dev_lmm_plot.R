p <- study_parameters(n1 = 11,
                      n2 = 2,
                      icc_pre_subject = 0.5,
                      effect_size = -1)


plot(p)

paras <- p
paras[is.na(paras)] <- 0

time <- get_time_vector(paras)

ES <- get_effect_size(paras)



y_c <- paras$fixed_intercept + paras$fixed_slope * time
y_tx <-  paras$fixed_intercept + (paras$fixed_slope + ) * time


y_c <- get_y(fixed_intercept = paras$fixed_intercept,
      fixed_slope = paras$fixed_slope,
      sigma_subject_intercept = paras$sigma_subject_intercept,
      sigma_cluster_intercept = paras$sigma_cluster_intercept,
      sigma_subject_slope = paras$sigma_subject_slope,
      sigma_cluster_slope = paras$sigma_cluster_slope,
      time = time,
      treatment = 0)

y_tx <-
    get_y(fixed_intercept = paras$fixed_intercept,
          fixed_slope = paras$fixed_slope + get_slope_diff(paras)/paras$T_end,
          sigma_subject_intercept = paras$sigma_subject_intercept,
          sigma_cluster_intercept = paras$sigma_cluster_intercept,
          sigma_subject_slope = paras$sigma_subject_slope,
          sigma_cluster_slope = paras$sigma_cluster_slope,
          time = time,
          treatment = 0)

d <- data.frame(y = c(y, y1),
                time = rep(time, 2),
                treatment = rep(c(0, 1), each = length(y)))
d$treatment <- factor(d$treatment, labels = c("Control", "Treatment"))

