# variance partitioning


#' Calculate the variance partitioning coefficient
#'
#' @param object An object created by \code{\link{study_parameters}}
#'
#' @details For partially nested studies, the VPC is calculated for the treatment group.
#'
#' @return a \code{data.frame} with class \code{plcp_VPC} containing the
#'  percentage of variance per level and time point. The column
#'  \code{between_clusters} is also the intraclass correlation for level three,
#'  i.e. the correlation between two subjects belonging to the same cluster at
#'  a specific time point. With random slopes in the model the variances per time point
#'  will be a quadratic function of time. \code{tot_var} is the
#'  percentage increase or decrease in total variance relative to baseline variance.
#'
#' The \code{plot} method returns a \code{ggplot2::ggplot} object.
#' @seealso \code{\link{plot.plcp_VPC}}
#'
#' @references Goldstein, H., Browne, W., & Rasbash, J. (2002).
#' Partitioning variation in multilevel models.
#' \emph{Understanding Statistics: Statistical Issues in Psychology, Education,
#' and the Social Sciences, 1}(4), 223-231.
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = 0.05,
#'                           var_ratio = 0.03)
#'
#' res <- get_VPC(paras)
#' res
#'
#' # Plot
#' plot(res)
#' @export
get_VPC <- function(object) {
     UseMethod("get_VPC")
}

#' @rdname get_VPC
#' @export
get_VPC.plcp <- function(object) {
    paras <- NA_to_zero(object)
     u0 <- paras$sigma_subject_intercept
     u1 <- paras$sigma_subject_slope
     v0 <- paras$sigma_cluster_intercept
     v1 <- paras$sigma_cluster_slope
     v01 <- v0 * v1 * paras$cor_cluster
     error <- paras$sigma_error

     u01 <- paras$cor_subject * u0 * u1

     time <- get_time_vector(paras)

     tot_var <- (u0^2 + 2*u01*time + u1^2*time^2 +
                     v0^2 + 2*v01*time + v1^2*time^2 + error^2)
     lvl3 <- (v0^2 + 2*v01*time + v1^2*time^2)/tot_var
     lvl2 <- (u0^2 + 2*u01*time + u1^2*time^2)/tot_var
     lvl1 <- error^2/tot_var

     tot_var <- tot_var/tot_var[1]
     res <- data.frame(time,
                       between_clusters = lvl3*100,
                       between_subjects = lvl2*100,
                       within_subjects = lvl1*100,
                       tot_var = (tot_var-1)*100)

     class(res) <- append("plcp_VPC", class(res))

     res
}
#' @export
get_VPC.plcp_multi <- function(object) {
    warning("Multiple study designs used, only the first is shown")
    get_VPC.plcp(object[1, ])
}


#' Plot method for \code{get_VPC}-objects
#'
#' @param x An object created with \code{\link{get_VPC}}
#' @param ... Optional arguments, currently ignored.
#'
#' @export

plot.plcp_VPC <- function(x, ...) {
     check_installed("ggplot2")
     res <- x
     res$tot_var <- NULL

     res <- stats::reshape(res, direction = "long",
                    varying = 2:4,
                    times = c("between_clusters",
                              "between_subjects",
                              "within_subjects"),
                    v.names = "proportion", timevar = "level")
     res$level <- factor(res$level, levels = c("between_clusters",
                                               "between_subjects",
                                               "within_subjects"),
                         labels = c("between-clusters (L3)",
                                    "between-subjects (L2)",
                                    "within-subjects (L1)"))

    p <- ggplot2::ggplot(res, ggplot2::aes_string("time", "proportion", color = "level", fill = "level")) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(title = "Variance partitioning",
               x = "Time point",
               y = "Proportion of total variance")

    if(requireNamespace("ggsci", quietly = TRUE)) {
       p <- p + ggsci::scale_fill_d3() +
            ggsci::scale_color_d3()
    }
    p
}

#' Print method for \code{get_vpc}-objects
#' @param x Object created with \code{link{get_VPC}}
#' @param digits Number of digits to print
#' @param ... Optional arguments
#' @method print plcp_VPC
#' @export
print.plcp_VPC <- function(x, digits = 2, ...) {
    cat("# Percentage (%) of total variance at each level and time point\n")
    print.data.frame(x, digits = digits, scientific = FALSE, ...)
    invisible(x)
}



# Standard deviations -----------------------------------------------------

#' Calculate the model implied standard deviations per time point
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param treatment \code{character}; either \code{"treatment"} or \code{"control"}.
#' Indicates for which group SDs should be calculated for. This only makes a difference
#' for 3-level partially nested designs.
#' @param n Optional; selects row n if \code{object} is a \code{data.frame} of
#' parameters
#'
#' @seealso \code{\link{get_VPC}}, \code{\link{get_correlation_matrix}}
#' @return \code{data.frame} with class \code{plcp_sds} containing the model
#' implied standard deviations per time point.
#'
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 6,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = 0.05,
#'                           var_ratio = 0.03)
#'
#' get_sds(paras)
#'
#' # plot
#' plot(get_sds(paras))
#'
get_sds <- function(object, treatment = "treatment", n = 1) {
    if(!treatment %in% c("treatment", "control")) stop("Wrong 'treatment', allowed options are: 'treatment' or 'control'", call. = FALSE)
     UseMethod("get_sds")
}

#' @export
get_sds.plcp <- function(object, treatment = "treatment", n = NULL) {
    .p <- NA_to_zero(object)
    .p <- prepare_paras(.p)
    if(treatment == "treatment") {
        .p <- .p$treatment
    } else if(treatment == "control") {
        .p <- .p$control
    }
    .p$retention <- NULL
    .p$n2 <- NULL
    .p <- .p[c("sigma_subject_intercept",
                 "cor_subject",
                 "cor_cluster",
                 "sigma_subject_slope",
                 "sigma_cluster_slope",
                 "sigma_cluster_intercept",
                 "sigma_error",
                 "n1",
                 "T_end")]

    res <- do.call(get_sds_, .p)

     class(res) <- append(c("plcp_sds"), class(res))

     res
}

#' @export
get_sds.plcp_multi <- function(object, treatment = "treatment", n = 1) {
    get_sds.plcp(as.plcp(object[n, ]), treatment = treatment)

}

get_sds_ <- function(sigma_subject_intercept,
                     cor_subject,
                     sigma_subject_slope,
                     sigma_cluster_intercept,
                     sigma_cluster_slope,
                     cor_cluster,
                     sigma_error,
                     n1,
                     T_end) {
     time <-  seq(0, T_end, length.out = n1)
     u0 <- sigma_subject_intercept
     u1 <- sigma_subject_slope
     u01 <- u0 * u1 * cor_subject
     v0 <- sigma_cluster_intercept
     v1 <- sigma_cluster_slope
     v01 <- v0 * v1 * cor_cluster

     error <- sigma_error
     sds <- sqrt((u0^2 + 2*u01*time + u1^2*time^2 + v0^2 +
                      2 * v01*time + v1^2*time^2 + error^2))

     sds_lvl2 <- sqrt((u0^2 + 2*u01*time + u1^2*time^2 + v0^2 + error^2))

     res <- data.frame(time = time,
                       SD_with_random_slopes = sds,
                       SD_no_cluster_random_slope = sds_lvl2,
                       SD_no_random_slopes = sqrt(u0^2 + v0^2 + error^2))

     res


}


#' Plot method for \code{get_sds}-objects
#' @param x An object of class \code{plcp_sds}.
#' @param ... Optional arguments.
#' @export
plot.plcp_sds <- function(x, ...) {
     check_installed("ggplot2")
    .res <- x
     cs <- .res$SD_no_random_slopes[1]

     res <- .res
     res$time <- round(res$time,1)

     p <- ggplot2::ggplot(res, ggplot2::aes_string("time", "SD_with_random_slopes")) +
         ggplot2::geom_hline(ggplot2::aes_string(color = "'Random slopes = 0'",
                         yintercept = "SD_no_random_slopes")) +
         ggplot2::geom_line(ggplot2::aes(color = "With random slopes")) +
         ggplot2::geom_point(ggplot2::aes(color = "With random slopes")) +
         ggplot2::scale_x_continuous(breaks = unique(res$time)) +
         ggplot2::labs(y = "SD", x = "Time point",
               title = "SD per time point",
               color = "Model")
     # facet_grid(~cor_cluster + cor_subject, labeller = label_both)

     #  if(!is.nulls(facet)) p + facet_wrap(facet)

     p

}

#' Print method for \code{get_sds}-objects
#' @param x An object of class \code{plcp_sds}.
#' @param ... Optional arguments.
#' @export
#' @method print plcp_sds
print.plcp_sds <- function(x, ...) {
    print.data.frame(x, digits = 2, ...)
}


#
# CORRELATIONS LVL 2
#

#' Calculate the subject-level (ICC) correlations among time points
#'
#' @param object An object created by \code{\link{study_parameters}}
#'
#' @return A \code{n1} x \code{n1} \code{matrix} with the marginal subject-level
#' correlations between time points.
#' @details The correlation between time point \eqn{T_i} and \eqn{T_{i+1}} within
#' the same subject is also called the intraclass correlation (ICC) at level two.
#' If the random slopes are non-zero this ICC change over time.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = 0.05,
#'                           var_ratio = 0.03)
#' get_correlation_matrix(paras)
get_correlation_matrix <- function(object) {
    UseMethod("get_correlation_matrix")
}

#' @export
get_correlation_matrix.plcp <- function(object) {
    paras <- NA_to_zero(object)


    u0 <- paras$sigma_subject_intercept
    u1 <- paras$sigma_subject_slope
    v0 <- paras$sigma_cluster_intercept
    v1 <- paras$sigma_cluster_slope
    v01 <- v0 * v1 * paras$cor_cluster
    error <- paras$sigma_error
    u01 <- paras$cor_subject * u0 * u1
    time <- get_time_vector(paras)

    n1 <- paras$n1
    n2 <- paras$n2
    sx2 <- sum( (time - mean(time))^2)/n1

    X <- matrix(c(rep(1, n1), time), ncol = 2)
    Z <- X
    D <- matrix(c(u0^2, u01,
                  u01, u1^2), ncol = 2)
    D2 <- matrix(c(v0^2, v01,
                   v01, v1^2), ncol = 2)
    V <- Z %*% D %*% t(Z) + Z %*% D2 %*% t(Z) + error^2*diag(n1)
    V <- cov2cor(V)

    time_rounded <- round(time, 1)
    dimnames(V) <- list(time_rounded, time_rounded)

    class(V) <- append(class(V), "plcp_ICC2")

    V
}


#' Plot method for \code{get_correlation_matrix}-objects
#'
#' @param x An object created with \code{\link{get_correlation_matrix}}
#' @param ... Optional arguments, currently ignored.
#'
#' @export


plot.plcp_ICC2 <- function(x, ...) {
    check_installed("ggplot2")
    res <- as.data.frame(x)

    breaks <- 1:ncol(res)
    res <- reshape(res, varying = breaks,
                   v.names = "cor",
                   idvar = "time1",
                   timevar = "time2",
                   direction = "long")
    res$time1 <- res$time1
    res$time2 <- res$time2
    res$cor2 <- round(res$cor, 2)

    break_labs <- as.numeric(dimnames(x)[[1]])

    p <- ggplot2::ggplot(res, ggplot2::aes_string("time1", "time2", color = "cor", fill = "cor")) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes_string(label = "cor2"), hjust = "center", color = "black") +
        ggplot2::scale_x_continuous(breaks = breaks, labels = break_labs) +
        ggplot2::scale_y_continuous(breaks = breaks, labels = break_labs) +
        ggplot2::labs(color = "Correlation", fill = "Correlation",
             x = "Time", y = "Time",
             title = "Subject-level correlation matrix") +
        ggplot2::theme_minimal()

    if(requireNamespace("viridis", quietly = TRUE)) {
        p <- p + viridis::scale_fill_viridis() +
            viridis::scale_color_viridis()
    }
    p

}

#' Print method for \code{get_correlation_matrix}-objects
#'
#' @param x An object created by \code{\link{get_correlation_matrix}}
#'
#' @param ... Optional arguments
#'
#' @method print plcp_ICC2
#' @export
print.plcp_ICC2 <- function(x, ...) {
    x <- unclass(x)
    print(round(x, 2), ...)
}

#' @rdname get_correlation_matrix
#' @export
get_correlation_matrix.plcp_multi <- function(object) {
    warning("Multiple study designs used, only the first is shown")
    get_correlation_matrix.plcp(object[1, ])
}


