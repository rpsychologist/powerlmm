# simulate_3lvl_data <- function(paras) {
#      UseMethod(".simulate_3lvl_data")
# }

# simulate_3lvl_data.list <- function(paras) {
#      do.call(.simulate_3lvl_data, paras)
# }
# simulate_3lvl_data.data.frame <- function(paras) {
#
#      lapply(1:nrow(paras), function(i) do.call(.simulate_3lvl_data, paras[i,]))
# }

create_cluster_index <- function(n2, n3) {
    if(length(n2) == 1) {
        cluster <- rep(1:n3, each = n2)
    } else {
        # if(is.null(n2_func)) {
        cluster <- lapply(seq_along(n2), function(i) rep((1:n3)[i], each = n2[i]))
        cluster <- unlist(cluster) # index clusters
        #} else if(n2_func == "runif") {
        #     n2 <- floor(runif(n3, n2[1], n2[2]))
        #}

    }
    cluster

}
.simulate_3lvl_data <- function (n1,
                           n2,
                           n3,
                           T_end,
                           fixed_intercept,
                           fixed_slope,
                           sigma_subject_intercept,
                           sigma_subject_slope,
                           sigma_cluster_intercept,
                           sigma_cluster_slope,
                           sigma_error,
                           shape,
                           cor_subject = 0,
                           cor_cluster = 0,
                           cor_within = 0,
                           dropout = NULL,
                           deterministic_dropout = NULL,
                           family = "gaussian",
                           ...) {

     # errors
     #if(!"MASS" %in% installed.packages()) stop("Package 'MASS' is not installed")

     # unbalanced n2
     #n2 <- unlist(n2)
   #  if(length(n2) != n3) stop("n2 and n3 do not mach")

     time <- seq(0, T_end, length.out = n1) # n1 measurements during the year

     n2_func <- names(n2)
     n2 <- unlist(n2)
     cluster <- create_cluster_index(n2, n3)

     subject <- rep(1:length(cluster), each = n1) # subject IDs
     tot_n2 <- length(cluster)

     # level-2 variance matrix
     Sigma_subject <- c(
          sigma_subject_intercept^2 ,
          sigma_subject_intercept * sigma_subject_slope * cor_subject,
          sigma_subject_intercept * sigma_subject_slope * cor_subject,
          sigma_subject_slope^2
     )
     Sigma_subject <- matrix(Sigma_subject, 2, 2) # variances

     # level-3 variance matrix
     Sigma_cluster <- c(
          sigma_cluster_intercept^2,
          sigma_cluster_intercept * sigma_cluster_slope * cor_cluster,
          sigma_cluster_intercept * sigma_cluster_slope * cor_cluster,
          sigma_cluster_slope^2
     )
     Sigma_cluster <- matrix(Sigma_cluster, 2, 2)

     # level 3-model
     cluster_lvl <-
          MASS::mvrnorm(sum(n3),
                  mu = c(0, 0),
                  Sigma = Sigma_cluster)

     if (is.null(dim(cluster_lvl))) {
          # if theres only one therapist
          cluster_b0 <- cluster_lvl[1]
          cluster_b1 <- cluster_lvl[2]
     } else {
          cluster_b0 <- cluster_lvl[, 1]
          cluster_b1 <- cluster_lvl[, 2]
     }

     cluster_b0 <- cluster_b0[cluster]
     cluster_b1 <- cluster_b1[cluster]

     # level 2-model
     subject_lvl <- MASS::mvrnorm(tot_n2, c(0, 0), Sigma_subject)

     b0 <- fixed_intercept + subject_lvl[, 1] + cluster_b0
     b1 <- fixed_slope + subject_lvl[, 2] + cluster_b1

     # linear predictor
     eta <-  b0[subject] + b1[subject] * time

     # Sample Y
     if(family == "gaussian") {

         # level-1 model
         sigma.y <- diag(n1)
         sigma.y <-
             sigma_error ^ 2 * cor_within ^ abs(row(sigma.y) - col(sigma.y)) # AR(1)

         # gen level-1 error
         error_sigma.y <- MASS::mvrnorm(tot_n2, rep(0, n1), sigma.y)

         y <- eta + c(t(error_sigma.y))
     } else if(family == "binomial") {
         y <- rbinom(tot_n2 * n1, 1, prob = plogis(eta))
     } else if(family == "poisson") {
         y <- rpois(tot_n2 * n1, lambda = exp(eta))
     }
     else if(family == "gamma") {
         y <- rgamma(tot_n2 * n1,
                     shape = shape,
                     rate = shape/exp(eta))
     } else if(family == "lognormal") {
         y <- rlnorm(tot_n2 * n1,
                     meanlog = eta,
                     sdlog = sigma_error)
     }

     df <-
          data.frame (y,
                      y_c = y,
                      time,
                      subject,
                      cluster = rep(cluster, each = n1),
                      fixed_intercept = fixed_intercept,
                      fixed_slope = fixed_slope,
                      subject_intercept = subject_lvl[, 1][subject],
                      subject_slope = subject_lvl[, 2][subject],
                      cluster_intercept = cluster_b0[subject],
                      cluster_slope = cluster_b1[subject],
                      miss = 0)

     df
}



