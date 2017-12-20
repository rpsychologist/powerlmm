# also add cohend = 0

nsim <- 10000
cores <- 31
cl <- parallel::makeCluster(cores)
# few_clusters ------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.1,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3,
                                                rate = 1/2),
                      partially_nested = FALSE,
                      cohend = -0.9)

cat("Sim 1\n")
res1 <- simulate(p, nsim = nsim, satterthwaite = TRUE, cores = cores, cl = cl)


## partially-nested
p <- update(p, partially_nested = TRUE, dropout = 0, n3 = 5, cohend = 0.6)
cat("Sim 2\n")
res2 <- simulate(p, nsim = nsim, satterthwaite = TRUE, cores = cores, cl = cl)



# more clusters ------------------------------------------------------------

p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 12,
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.1,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3,
                                                rate = 1/2),
                      partially_nested = FALSE,
                      cohend = 0.5)
cat("Sim 3\n")
res3 <- simulate(p, nsim = nsim, satterthwaite = TRUE, cores = cores, cl = cl)

# unequal clusters
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(5,5,5,6,7,10,10,15,20,20,25,30),
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.1,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3,
                                                rate = 1/2),
                      partially_nested = FALSE,
                      cohend = 0.5)
cat("Sim 4\n")
res4 <- simulate(p, nsim = nsim, satterthwaite = TRUE, cores = cores, cl = cl)

# Random clusters
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(func = rnorm(18, 5, 15), replace = 4),
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.1,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3,
                                                rate = 1/2),
                      partially_nested = FALSE,
                      cohend = 0.7)

cat("Sim 5\n")
res5 <- simulate(p, nsim = nsim, satterthwaite = TRUE, cores = cores, cl = cl)
parallel::stopCluster(cl)
saveRDS(list(res1, res2, res3, res4, res5), file = "simres.rds")
