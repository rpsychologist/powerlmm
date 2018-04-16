## ICC SLOPE

#' Calculate the amount of slope variance at the third level
#'
#' @param object An object created by \code{\link{study_parameters}}.
#' @param ... Optional named arguments.
#' @return Returns the proportion of slope variance at the third level as
#' a numeric vector. \code{NA} is returned for models with no slope variance
#' as either level two or three.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           sigma_subject_intercept = 1.2,
#'                           sigma_subject_slope = 0.2,
#'                           sigma_cluster_intercept = 0,
#'                           sigma_cluster_slope = 0.2,
#'                           sigma_error = 1.2,
#'                           cohend = -0.8)
#'
#' get_ICC_slope(paras)
get_ICC_slope <- function(object, ...) {
    UseMethod("get_ICC_slope")
}

#' @export
#'
get_ICC_slope.default <- function(object = NULL, u1, v1, ...) {
    u1[is.na(u1)] <- 0
    x <- v1^2/(v1^2 + u1^2)
    x[v1 == 0 & u1 == 0] <- 0
    x
}

#' @export
get_ICC_slope.plcp <- function(object, ...) {
    get_ICC_slope.default(u1 = object$sigma_subject_slope,
                      v1 = object$sigma_cluster_slope)
}

#' @export
get_ICC_slope.plcp_multi <- function(object, ...) {
    get_ICC_slope.plcp(object)
}


#' Calculates the ratio of the slope variance to the within-subjects error variance
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param ... Optional arguments.
#' @return Returns the ratio of the total slope variance to the within-subject error
#'  as a numeric vector.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           sigma_subject_intercept = 1.2,
#'                           sigma_subject_slope = 0.2,
#'                           sigma_cluster_intercept = 0,
#'                           sigma_cluster_slope = 0.2,
#'                           sigma_error = 1.2,
#'                           cohend = -0.8)
#'
#' get_var_ratio(paras)
get_var_ratio <- function(object, ...) {
    UseMethod("get_var_ratio")
}
#' @export
get_var_ratio.default <- function(object=NULL, v1, u1, error, ...) {
     #if(!(is.na(v1) && is.na(u1))) {
         v1[is.na(v1) & !is.na(u1)] <- 0
         u1[is.na(u1) & !is.na(v1)] <- 0
     #}

     (v1^2 + u1^2)/(error^2)
}
#' @export
get_var_ratio.plcp <- function(object, ...) {
    get_var_ratio.default(v1 = object$sigma_cluster_slope,
                            u1 = object$sigma_subject_slope,
                            error = object$sigma_error)
}
#' @export
get_var_ratio.plcp_multi <- function(object, ...) {
    get_var_ratio.plcp(object)
}

##
## ICC_pre_subjects
##


#' Calculate the subject-level ICC at pretest
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param ... Optional named arguments.
#'
#' @return Returns the proportion of baseline variance at the subject level (which also includes cluster-level variance),
#' as a numeric vector.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           sigma_subject_intercept = 1.2,
#'                           sigma_subject_slope = 0.2,
#'                           sigma_cluster_intercept = 0.5,
#'                           sigma_cluster_slope = 0.2,
#'                           sigma_error = 1.2,
#'                           cohend = -0.8)
#'
#' get_ICC_pre_subjects(paras)
get_ICC_pre_subjects <- function(object, ...) {
    UseMethod("get_ICC_pre_subjects")
}
#' @export
get_ICC_pre_subjects.default <- function(object = NULL, u0, v0, error, ...) {
    v0 <- ifelse(is.na(v0), 0, v0)

    (u0^2 + v0^2)/(u0^2 + v0^2 + error^2)
}
#' @export
get_ICC_pre_subjects.plcp <- function(object, ...) {
    get_ICC_pre_subjects.default(u0 = object$sigma_subject_intercept,
                            v0 = object$sigma_cluster_intercept,
                            error = object$sigma_error)
}
#' @export
get_ICC_pre_subjects.plcp_multi <- function(object, ...) {
    get_ICC_pre_subjects.plcp(object)
}

##
## ICC_pre_clusters
##

#' Calculate the amount of baseline variance at the cluster level
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param ... Optional named arguments.
#' @details
#' The proportion of variance at the cluster level at baseline can
#' be interpreted as the correlation between two subjects belonging to the same
#' cluster.
#'
#' @return Returns the proportion of baseline variance at the cluster level,
#' as a numeric vector.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           sigma_subject_intercept = 1.2,
#'                           sigma_subject_slope = 0.2,
#'                           sigma_cluster_intercept = 0.5,
#'                           sigma_cluster_slope = 0.2,
#'                           sigma_error = 1.2,
#'                           cohend = -0.8)
#'
#' get_ICC_pre_clusters(paras)
get_ICC_pre_clusters <- function(object, ...) {
    UseMethod("get_ICC_pre_clusters")
}
#' @export
get_ICC_pre_clusters.default <- function(object = NULL, u0, v0, error, ...) {
    v0^2/(u0^2 + v0^2 + error^2)
}
#' @export
get_ICC_pre_clusters.plcp <- function(object, ...) {
    get_ICC_pre_clusters.default(u0 = object$sigma_subject_intercept,
                            v0 = object$sigma_cluster_intercept,
                            error = object$sigma_error)
}
#' @export
get_ICC_pre_clusters.plcp_multi <- function(object, ...) {
    get_ICC_pre_clusters.plcp(object)
}



var_T <- function(n1, T_end) {
    time <- seq(0, T_end, length.out = n1)
    sx <- sum( (time - mean(time))^2)/n1
}
get_ICC_aov <- function(u0, u1, u01 = 0, t, v1, v0 = 0, v01 = 0, error) {
    UseMethod("get_ICC_aov")
}
get_ICC_aov.default <- function(u0, u1, u01 = 0, t, v1, v0 = 0, v01 = 0, error) {

    f2 <- (u0^2 + 2*u01*t + t^2*u1^2 + v0^2 + 2*v01*t + t^2 * v1^2 + error^2)
    ICC <- (v0^2 + 2*v01*t + t^2 * v1^2)/f2

    ICC
}
get_ICC_aov.list <- function(x) {
    args <- list()
    args$u0 <- x$sigma_subjects_intercept
    args$u1 <- x$sigma_subject_slope
    args$u01 <- x$cor_subject * args$u0 * args$u1
    args$t <- x$T_end
    args$v1 <- x$sigma_cluster_slope
    args$v0 <- x$sigma_cluster_intercept
    args$v01 <- x$cor_cluster *  args$v0 * args$v1
    args$error <- x$sigma_error

    do.call(get_ICC_aov.default, args)
}

get_time_vector <- function(paras) {
     seq(0, paras$T_end, length.out = paras$n1)
}


# logic -------------------------------------------------------------------

NA_or_zero <- function(x) {
    x == 0 | is.na(x)
}

# convert NA paras to 0
# used e.g. in get_SDS()
NA_to_zero <- function(object) {
    ind <- c("sigma_subject_intercept", "sigma_subject_slope", "cor_subject",
             "sigma_cluster_intercept", "sigma_cluster_slope", "cor_cluster",
             "sigma_error")

    for(i in ind) {
        x <- object[[i]]
        object[[i]][is.na(x)] <- 0
    }

    object
}

# elapsed time --------------------------------------------------------------------
elapsed_time <- function(object, ...) {
     UseMethod("elapsed_time")
}

elapsed_time.plcp_multi_sim <- function(object) {
     x <- lapply(object, function(d) {
          n1 <- d$paras$n1
          n2 <- d$paras$n2
          n3 <- d$paras$n3
          time <- d$time

          data.frame("n1" = n1, "n2" = n2, "n3" = n3, "elapsed_time" = time)

     })
     x <- do.call(rbind, x)


     class(x) <- append("elapsed_time", class(x))
     attr(x, "nsim") <- object[[1]]$nsim
     x
}
print.elapsed_time <- function(x, ...) {
     print.data.frame(x)
     cat("# nsim: ", attr(x, "nsim"))
}


# check if approx equal --------------------------------------------------------------------
is_approx <- function(x, y) {
    abs(x - y) < .Machine$double.eps^0.5
}


# save multi-simulations --------------------------------------------------------
save_res <- function(res, i) {
     output_dir <- format(Sys.time(), "%Y%m%d_%H%M")
     output_dir <- paste("save/", output_dir, sep = "")
     dir.create(output_dir)
     f <- paste(output_dir, "/sim", i, ".rds", sep = "")

}


# load saved multi-sims ---------------------------------------------------------
load_sims <- function(folder) {
     files <- dir(folder)
     files <- paste(folder, files, sep ="")
     x <- lapply(files, readRDS)
     class(x) <- append(class(x), "multi_sim")

     x
}


# proportion mising per time point and group
sum_missing_tx_time <- function(.d) {
    d <- with(.d, tapply(y, list(treatment, time), function(x) mean(is.na(x))))
    d_c <- data.frame(treatment = 0,
                      time = unique(.d$time),
                      missing = d[1,])
    d_tx <- data.frame(treatment = 1,
                       time = unique(.d$time),
                       missing = d[2,])

    rbind(d_c, d_tx)

}


# Plot design
#' Plot method for \code{study_parameters}-objects
#' @param x An object of class \code{plcp}.
#' @param n specifies which row \code{n} should be used if \code{object}
#'  is a \code{data.frame} containing multiple setups.
#' @param type indicated what plot to show. If \code{effect} the plot showing the treatment groups
#' change over time will be shown, if \code{dropout} the missing data pattern will be shown,
#' if \code{both} both plots will be shown.
#'
#' @param ... Optional arguments.
#' @export
plot.plcp <- function(x, n = 1, type = "both", ...) {
    check_installed("ggplot2")
    paras <- x
     if(is.data.frame(paras)) {
          paras <- as.list(paras[n,])
          paras <- do.call(study_parameters, paras)
          #class(paras) <- append(c("plcp"), class(paras))
     }

     time <- get_time_vector(paras)


     y <- paras$fixed_intercept + paras$fixed_slope * time
     y1 <-  paras$fixed_intercept + (paras$fixed_slope + get_slope_diff(paras)/paras$T_end) * time

     d <- data.frame(y = c(y, y1),
                     time = rep(time, 2),
                     treatment = rep(c(0, 1), each = length(y)))
     d$treatment <- factor(d$treatment, labels = c("Control", "Treatment"))

     p1 <- ggplot2::ggplot(d, ggplot2::aes_string("time", "y", color = "treatment")) +
         ggplot2::geom_line(show.legend = TRUE) +
         ggplot2::geom_point(show.legend = FALSE) +
         ggplot2::labs(title = "Treatment effects", y = "Outcome", x = "Time point",
               subtitle = paste("Difference at endpoint is equal to Cohen's d =", paras$cohend),
               caption = "N.B Cohen's d is calculated using baseline standard deviations")

     d <- simulate_data(paras)
     d <- sum_missing_tx_time(d)
     d$treatment <- factor(d$treatment, labels = c("Control", "Treatment"))



     #theoretical_missing <- get_dropout(update(paras, n1 = 100))
     theoretical_missing <- get_dropout(paras)

     theoretical_missing_tx <- theoretical_missing[ , c("time", "treatment")]
     colnames(theoretical_missing_tx)[2] <- "missing"
     theoretical_missing_tx$treatment <- 1
     theoretical_missing_cc <- theoretical_missing[ , c("time", "control")]
     theoretical_missing_cc$treatment <- 0
     colnames(theoretical_missing_cc)[2] <- "missing"
     theoretical_missing <- rbind(theoretical_missing_cc,
                                  theoretical_missing_tx)
     theoretical_missing$treatment <- factor(theoretical_missing$treatment,
                                             labels = c("Control", "Treatment"))


    p2 <- ggplot2::ggplot(d, ggplot2::aes_string("time", "missing", color = "treatment", group = "treatment")) +
            ggplot2::geom_point() +
            ggplot2::geom_line(data = theoretical_missing,
                        ggplot2::aes_string("time", "missing",
                          color = "treatment", group = "treatment"),
                      linetype = "dashed") +
        ggplot2::labs(title = "Dropout", y = "Proportion dropout", x = "Time point") +
        ggplot2::ylim(0,1)

     if(type == "both") {
         check_installed("gridExtra")
        return(gridExtra::grid.arrange(p1, p2, ncol=2))
     } else if(type == "effect") {
         return(p1)
     } else if(type == "dropout") {
         return(p2)
     }

}

#' @export
plot.plcp_multi <- function(x, n = 1, type = "both", ...) {
    plot.plcp(x, n = n, type = type)
}

# power curve


#' Create a power table for a combination of parameter values
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param n2 A vector of n2 values
#' @param ... Optional named arguments. Up to two extra arguments can be compared.
#' When used together with the plot method, the first argument will be grouped by
#' color and the second by facets.
#' @param df Either "between" or "satterth" for Satterthwaite's DF approximation.
#' Also accepts a \code{numeric} value which will be used as DF. See \code{\link{get_power}}
#' @param alpha The alpha level, defaults to 0.05.
#' @param R An \code{integer} indicating how many realizations to base power on.
#' Useful when dropout or cluster sizes are sampled (i.e. are random variables).
#' @param cores An \code{integer} indicating how many CPU cores to use.
#'
#'
#' @return A \code{data.frame} with class \code{plcp_power_table}.
#' @export
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 6,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           var_ratio = 0.03,
#'                           icc_slope = 0.05,
#'                           cohend = -0.8)
#'
#' # increase only n2
#' x <- get_power_table(paras, n2 = 10:15)
#' plot(x)
#'
#' # Compare two parameters
#' x <- get_power_table(paras, n2 = 10:15, n3 = 6:8)
#' plot(x)
#'
#' # Compare impact of three parameters
#' x <- get_power_table(paras, n2 = seq(3, 25, by = 3),
#'                             n3 = c(3,6,9),
#'                             icc_slope = c(0, 0.05, 0.1))
#' plot(x)
get_power_table <- function(object, n2, ..., df = "between", alpha = 0.05, R = 1L, cores = 1L) {

    paras <- object
    arg <- list(...)
    updateProgress <- arg$updateProgress
    arg$updateProgress <- NULL
    if(length(arg) > 2) warning("Too many extra arguments, only the first 2 will be used")
    arg <- arg[1:2]
    extra_args <- c(names(arg)[1], names(arg)[2])
    arg$object <- paras
    arg$n2 <- n2
    arg[which(names(arg) == "")] <- NULL
    arg[which(is.na(names(arg)))] <- NULL

    paras <- do.call(update.plcp, arg)
    if(inherits(paras, "plcp")) stop("Only 1 power value requested, add more parameter values.", call. = FALSE)

    res <- get_power(paras, updateProgress = updateProgress, df = df, alpha = alpha, R = R, cores = cores)

    tmp <- paras
    tmp$icc_slope <- get_ICC_slope(paras)
    tmp$icc_pre_cluster <- get_ICC_pre_clusters(paras)
    tmp$icc_pre_subject <- get_ICC_pre_subjects(paras)
    tmp$var_ratio <- get_var_ratio(paras)

    if(length(extra_args) == 0) {
        extra_args <- character(0)
    } else {
        extra_args <- extra_args[!is.na(extra_args) & extra_args != ""]
    }

    tmp <- tmp[, c("n2", extra_args), drop = FALSE]

    tot_n <- lapply(res$tot_n, colMeans)
    tot_n <- do.call(rbind, tot_n)
    tot_n <- as.data.frame(tot_n)
    res <- cbind(tmp, power = unlist(res$power),
                 tot_n = tot_n$total)
    res$dropout <- "with missing"

    if(is.list(paras$dropout)) {
        res2 <- get_power(update(paras, dropout = 0), df = df, alpha = alpha, R = R, cores = cores)
        tot_n <- lapply(res2$tot_n, colMeans)
        tot_n <- do.call(rbind, tot_n)
        tot_n <- as.data.frame(tot_n)
        res2 <- cbind(tmp, power = unlist(res2$power),
                      tot_n = tot_n$total)
        res2$dropout <- "no missing"

        res <- rbind(res, res2)
    } else {
        res$dropout <- "no missing"
    }
    #res$tot_n <- paras$n2  * paras$n3

    if("effect_size" %in% colnames(res)) {
        res$effect_size <- get_effect_size.plcp_multi(res)$ES
    }



    for(i in extra_args) {
        res[, i] <- factor(res[, i])
    }

    res$dropout <- factor(res$dropout)
    class(res) <- append("plcp_power_table", class(res))


    res

}
#' Plot method for \code{get_power_table}-objects
#' @param x An object of class \code{plcp_power_table}.
#' @param ... Optional arguments.
#' @export
plot.plcp_power_table <- function(x, ...) {
    check_installed("ggplot2")
    .d <- x
    x <- .d[, which(!colnames(.d) %in% c("n2","power", "dropout", "tot_n")), drop = FALSE]
     .d$id <- interaction(cbind(x, .d$dropout))


    facet <- colnames(x)[2]
    color <- colnames(x)[1]
    if(is.na(color)) color <- NULL
    .d$dropout <- factor(.d$dropout)
    p <- ggplot2::ggplot(.d, ggplot2::aes_string("tot_n", "power", group = "id", color = color, linetype = "dropout")) +
        ggplot2::geom_line() +
        ggplot2::geom_point(ggplot2::aes_string(alpha = "dropout")) +
        ggplot2::theme_minimal() +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(10)) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(10)) +
        ggplot2::scale_alpha_manual(values = c("with missing" = 0, "no missing" = 1), guide = FALSE) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = "dotted") +
        ggplot2::labs(linetype = "missing data",
             color = colnames(.d)[2],
             x = "total number of subjects in study",
             y = "Power",
             title = "Power curves")


    if(requireNamespace("ggsci", quietly = TRUE)) p <- p + ggsci::scale_color_d3()
    if(!is.na(facet)) p <- p + ggplot2::facet_grid(as.formula(paste("~", facet)), labeller = "label_both")

    p
}



## drop out
format_dropout <- function(dropout) {
    dropout <- dropout * 100
    dropout <- round(dropout, 2)
    dropout <- paste(format(dropout, nsmall = 0), collapse = ",")
    dropout <- paste(dropout, "(%)")
}
make_list_weibull <- function(x) {
    if(all(class(x) == c("function", "plcp_weibull"))) {
        x <- list(x)
    }

    x
}

# monte carlo error -------------------------------------------------------
#' Calculate the Monte Carlo standard error of the empirical power estimates
#'
#' Returns the expected simulation error for a study design. Indicates how many
#' simulation that are needed for a desired precision in the empirical power
#' estimates.
#'
#' @param object An object created by \code{\link{get_power}}
#' @param nsim A \code{numeric} indicating the number of simulations
#' @param power \emph{Optional}. A \code{numeric} indicating the empirical power.
#' @param ... Currently not used.
#' Used when \code{object} is \code{NULL}.
#'
#' @return A \code{data.frame} with the estimated power, expected standard error
#'  of the simulated power estimate, and the 95 \% CI of the estimate.
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
#'                           var_ratio = 0.03,
#'                           icc_slope = 0.05,
#'                           cohend = -0.8)
#'
#' x <- get_power(paras)
#' get_monte_carlo_se(x, nsim = 1000)
#'
#' # Without an object
#' get_monte_carlo_se(power = 0.8, nsim = 1000)
get_monte_carlo_se <- function(object, nsim, power, ...) {
    UseMethod("get_monte_carlo_se")
}
#' @export
get_monte_carlo_se.default <- function(object, nsim, power, ...) {
    get_monte_carlo_se_(p = power, nsim)
}

get_monte_carlo_se_ <- function(p, nsim) {
    se <- sqrt((p * (1-p)/nsim))

    res <- data.frame(power = p,
               se = se,
               lwr_95 = p - 1.96*se,
               upr_95 = p + 1.96*se)

    class(res) <- append("plcp_mc_se", class(res))

    res
}

# report MCSE based on Gaussian approx
get_monte_carlo_se_gaussian <- function(x) {
    sd(x)/sqrt(length(x))
}

#' @rdname get_monte_carlo_se
#' @export
get_monte_carlo_se.plcp_power_3lvl <- function(object, nsim, ...) {
    p <- object$power
    get_monte_carlo_se_(p, nsim)
}
#' @rdname get_monte_carlo_se
#' @export
get_monte_carlo_se.plcp_power_2lvl <- function(object, nsim, ...) {
    p <- object$power
    get_monte_carlo_se_(p, nsim)
}
#' Print method for \code{get_monte_carlo_se}-objects
#' @param x An object created with \code{\link{get_monte_carlo_se}}.
#' @param digits The number of digits to print.
#' @param ... Optional arguments.
#' @method print plcp_mc_se
#' @export
print.plcp_mc_se <- function(x, digits = 2, ...) {
    print.data.frame(x, digits = digits, ...)
}



# Shiny -------------------------------------------------------------------
#' Launch powerlmm's Shiny web application
#'
#' This Shiny application provides the basic functionality of the \pkg{powerlmm}-package
#' in a user friendly web application.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' library(shiny)
#' shiny_powerlmm()
#' }
shiny_powerlmm <- function() {
    check_installed(c("shiny", "shinydashboard", "viridis", "ggsci", "ggplot2", "tidyr"))
    shiny::runApp(system.file("shiny_powerlmm", package = "powerlmm"))
}

# check installed packages
check_installed <- function(x) {
   pkg <- lapply(x, function(y) requireNamespace(y, quietly = TRUE))
   pkg <- unlist(pkg)

   x <- x[!pkg]

   if(length(x) > 0) stop(paste("Please install packages: ", paste(x, collapse = ", ")))
}
