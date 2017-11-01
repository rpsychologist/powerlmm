library(powerlmm)
library(methods)


# setup -------------------------------------------------------------------
# two-level
p2 <- study_parameters(n1 = 11,
                       n2 = 10,
                       n3 = c(4, 8),
                       T_end = 10,
                       fixed_intercept = 37,
                       fixed_slope = -0.64,
                       sigma_subject_intercept = 2.8,
                       sigma_subject_slope = 0.4,
                       sigma_cluster_intercept = 0,
                       cor_subject = -0.4,
                       icc_slope = 0,
                       sigma_error = 2.6,
                       dropout = c(0, dropout_weibull(0.3, 0.5)),
                       cohend = c(0.5, 0.8))


# three-level, only random slope

n2 <- c(10, unequal_clusters(4,5,6,7,25))
p3s <- study_parameters(n1 = 11,
                        n2 = n2,
                        n3 = c(4, 8),
                        T_end = 10,
                        fixed_intercept = 37,
                        fixed_slope = -0.64,
                        sigma_subject_intercept = 2.8,
                        sigma_subject_slope = 0.4,
                        sigma_cluster_intercept = 0,
                        cor_subject = -0.4,
                        icc_slope = 0.05,
                        sigma_error = 2.6,
                        partially_nested = FALSE,
                        dropout = c(0, dropout_weibull(0.3, 0.5)),
                        cohend = c(0.5, 0.8))

p3spn <- study_parameters(n1 = 11,
                         n2 = 10,
                         n3 = c(4, 8),
                         T_end = 10,
                         fixed_intercept = 37,
                         fixed_slope = -0.64,
                         sigma_subject_intercept = 2.8,
                         sigma_subject_slope = 0.4,
                         sigma_cluster_intercept = 0,
                         cor_subject = -0.4,
                         icc_slope = 0.05,
                         sigma_error = 2.6,
                         partially_nested = TRUE,
                         dropout = c(0, dropout_weibull(0.3, 0.5)),
                         cohend = c(0.5, 0.8))


# three-level all params
p3a <- study_parameters(n1 = 11,
                        n2 = 10,
                        n3 = c(4, 8),
                        T_end = 10,
                        fixed_intercept = 37,
                        fixed_slope = -0.64,
                        sigma_subject_intercept = 2.8,
                        sigma_subject_slope = 0.4,
                        icc_pre_cluster = 0.1,
                        cor_subject = -0.4,
                        cor_cluster = -0.3,
                        icc_slope = 0.05,
                        sigma_error = 2.6,
                        partially_nested = FALSE,
                        dropout = c(0, dropout_weibull(0.3, 0.5)),
                        cohend = c(0.5, 0.8))

p3apn <- study_parameters(n1 = 11,
                          n2 = 10,
                          n3 = c(4, 8),
                          T_end = 10,
                          fixed_intercept = 37,
                          fixed_slope = -0.64,
                          sigma_subject_intercept = 2.8,
                          sigma_subject_slope = 0.4,
                          icc_pre_cluster = 0.1,
                          cor_subject = -0.4,
                          cor_cluster = -0.3,
                          icc_slope = 0.05,
                          sigma_error = 2.6,
                          partially_nested = TRUE,
                          dropout = c(0, dropout_weibull(0.3, 0.5)),
                          cohend = c(0.5, 0.8))



# simulate ----------------------------------------------------------------
nsim <- 5000
cores <- 6

cat("## p2 (1/5)")
res_p2 <- simulate(p2, nsim = nsim, formula = "y ~ time*treatment + (1 + time | subject)", cores = cores, satterthwaite = TRUE)

cat("## p3s (2/5)")
res_p3s <- simulate(p3s, nsim = nsim, formula = "y ~ time*treatment + (1 + time | subject) + (0 + time | cluster)", cores = cores, satterthwaite = TRUE)

cat("## p3spn (3/5)")
res_p3spn <- simulate(p3spn, nsim = nsim, formula = "y ~ time*treatment + (1 + time | subject) + (0 + treatment:time | cluster)", cores = cores, satterthwaite = TRUE)

cat("## p3a (4/5)")
res_p3a <- simulate(p3a, nsim = nsim, formula = "y ~ time*treatment + (1 + time | subject) + (1 + time | cluster)", cores = cores, satterthwaite = TRUE)

cat("## p3apn (5/5)")
res_p3apn <- simulate(p3apn, nsim = nsim, formula = "y ~ time*treatment + (1 + time | subject) + (0 + treatment + treatment:time | cluster)", cores = cores, satterthwaite = TRUE)

save(res_p2, res_p3s, res_p3spn, res_p3a, res_p3apn, file = "sims.Rdata")

