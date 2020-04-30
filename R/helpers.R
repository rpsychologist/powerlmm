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
get_ICC_slope.plcp_crossed <- function(object, interaction = FALSE) {
    u1 <- object$sigma_subject_slope
    v2 <- object$sigma_cluster_slope
    v3 <- object$sigma_cluster_slope_crossed

    if(interaction) {
        # the propoprtion of 3-lvl slope variance
        # caused by the interaction
        v3^2/(v2^2 + v3^2)
    } else {
        # proportion of total re slope variance
        # at the 3rd level
        (v2^2 + v3^2)/(u1^2 + v2^2 + v3^2)
    }
   
}

#' @export
get_ICC_slope.plcp_multi <- function(object, ...) {
    get_ICC_slope.plcp(object)
}
#' @export
get_ICC_slope.plcp_multi_crossed <- function(object, ...) {
    get_ICC_slope.plcp_crossed(object, ...)
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
get_ICC_pre_subjects.plcp_crossed <- function(object, ...) {
    u0 <- object$sigma_subject_intercept
    v0 <- object$sigma_cluster_intercept
    v1 <- object$sigma_cluster_intercept_crossed
    error <- object$sigma_error

    (u0^2 +  v0^2 + v1^2)/(u0^2 + v0^2 + v1^2 + error^2)
}

#' @export
get_ICC_pre_subjects.plcp_multi <- function(object, ...) {
    get_ICC_pre_subjects.plcp(object)
}
#' @export
get_ICC_pre_subjects.plcp_multi_crossed <- function(object, ...) {
    get_ICC_pre_subjects.plcp_crossed(object)
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
get_ICC_pre_clusters.plcp_crossed <- function(object, interaction = FALSE) {
    u0 <- object$sigma_subject_intercept
    v0 <- object$sigma_cluster_intercept
    v1 <- object$sigma_cluster_intercept_crossed
    error <- object$sigma_error

    if(interaction) {
        v1^2/(v0^2 + v1^2)
    } else {
        (v0^2 + v1^2)/(u0^2 + v0^2 + v1^2 + error^2)
    }
}
#' @export
get_ICC_pre_clusters.plcp_multi <- function(object, ...) {
    get_ICC_pre_clusters.plcp(object)
}
#' @export
get_ICC_pre_clusters.plcp_multi_crossed <- function(object, ...) {
    get_ICC_pre_clusters.plcp_crossed(object)
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

# Calculate RE intervals
#
# @param d simulate_data data.frame
# @param var variable to summarise
# @param treatment treatment group indicator
#
# @keywords internal
eta_sum_d <- function(d, var, treatment) {
    x <- lapply(unique(d$time), function(i) {
        x <- eta_sum(d[d$treatment == treatment & d$time == i, var])
        x <- as.data.frame(x)
        x$treatment <- treatment
        x$time <- i

        x
    })
    x <- do.call(rbind, x)

    x
}

# Reshape eta_sum output from wide to long
#
# @keywords internal
# @param x eta_sum_* data.frame
#
# @return data.frame in long format
reshape_eta_sum <- function(x) {

    tmp <- x[, !colnames(x) %in% c("var","mean", "sd", "Q50")]
    tmp <- lapply(list(c("Q0.5", "Q99.5"),
                       c("Q2.5", "Q97.5"),
                       c("Q10", "Q90"),
                       c("Q25", "Q75")),
                  function(Q) {
                      d <- tmp[, c("treatment","time", Q)]
                      data.frame(treatment = x$treatment,
                                 time = x$time,
                                 min = x[, Q[1]],
                                 max = x[, Q[2]],
                                 width = switch(Q[1],
                                                "Q0.5" = 0.99,
                                                "Q2.5" = 0.95,
                                                "Q10" = 0.8,
                                                "Q25" = 0.5)
                      )
                  })
    tmp <- do.call(rbind, tmp)
    tmp$width <- factor(tmp$width,
                        levels = rev(c(0.5, 0.8, 0.95, 0.99)))

    tmp
}

.rbind_lists <- function(args, func = NULL) {
    x <- lapply(seq_along(args), function(i) {
        if(is.null(func)) {
            tmp <- args[[i]]
        } else {
            tmp <- func(args[[i]])
        }

        tmp$var <- names(args)[i]

        tmp
    })
    x <- do.call(rbind, x)

    x
}
.get_facet_lims <- function(d, var_names, min_cols, max_cols, trim = c(0, 1)) {
    lims <- lapply(var_names, function(x) {
        tmp <- d[d$var == x, ]

        if(is.null(dim(tmp[, min_cols]))) {
            mean <- c(quantile(tmp[, min_cols], trim[1]),
                 quantile(tmp[, max_cols], trim[2]))
        } else {
            mean <- c(apply(tmp[, min_cols], 2, quantile, probs = trim[1]),
                     apply(tmp[, max_cols], 2, quantile, probs  = trim[2]))
        }

        data.frame(var = x,
                   #mean = c(min(tmp[, min_cols]), max(tmp[, max_cols])),
                   mean = mean,
                   treatment = "Treatment",
                   time = 0)
    })
    lims <- do.call(rbind, lims)

    lims
}


# Calc linera predictor level 2 and 3
#
# @param d data.frame for tx group
# @treatment treatment indicator
#
# @return a data.frame with the extra cols mu2 and mu3
# @keywords internal
.calc_mu <- function(d, p, treatment = 1) {
    d$treatment <- treatment

    if(treatment == 1) {
        slope_diff <- get_slope_diff(p)/p$T_end
    } else slope_diff <- 0

    # partialyl nested
    if(p$partially_nested & treatment == 0) {
        d$cluster_intercept <- 0
        d$cluster_slope <- 0
    }

    # level 2
    d$mu2 <- with(d, (p$fixed_intercept + subject_intercept + cluster_intercept) + (p$fixed_slope + slope_diff + subject_slope + cluster_slope) * time)

    # level 3
    d$mu3 <- with(d, p$fixed_intercept + cluster_intercept + (p$fixed_slope  + slope_diff + cluster_slope) * time)

    d
}

.mu_vec_to_long <- function(x, RE_level,
                            var1 = "mu1_vec",
                            var2 = "mu2_vec",
                            var3 = "mu3_vec",
                            tx_var = "y2",
                            level1_func = .sample_level1_nested,
                            ...) {
    res1 <- NULL
    res2 <- NULL
    res3 <- NULL

    if(any(RE_level == 1)) {
        x1 <- level1_func(x$paras,
                                    ...)
        res1 <- lapply(1:nrow(x1[["y"]]), function(i) {
            tmp <- data.frame(y = x1[[var1]][[i]])
            tmp$treatment <- x1[["y"]][i, "treatment"]
            tmp$time <- x1[["y"]][i, "time"]
            tmp$var <- "within-subject"
            tmp
        })
        res1 <- do.call(rbind, res1)
    }
    if(any(RE_level == 2)) {
        res2 <- lapply(1:nrow(x[[tx_var]]), function(i) {
            tmp <- data.frame(y = x[[var2]][[i]])
            tmp$treatment <- x[[tx_var]][i, "treatment"]
            tmp$time <- x[[tx_var]][i, "time"]
            tmp$var <- "subject"
            tmp
        })
        res2 <- do.call(rbind, res2)
    }
    if(any(RE_level == 3)) {
        res3 <- lapply(1:nrow(x[[tx_var]]), function(i) {
            tmp <- data.frame(y = x[[var3]][[i]])
            tmp$treatment <- x[[tx_var]][i, "treatment"]
            tmp$time <- x[[tx_var]][i, "time"]
            tmp$var <- "cluster"
            tmp
        })
        res3 <- do.call(rbind, res3)
    }

    res <- rbind(res1, res2, res3)
    res$var <- factor(res$var, labels = c("Within-subject","Subject", "Cluster"),
                      levels = c("within-subject","subject", "cluster"))
    res$treatment <- factor(res$treatment, labels = c("Control", "Treatment"))

    res
}

.plot_dropout <- function(paras) {
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

   ggplot2::ggplot(d, ggplot2::aes_string("time", "missing", color = "treatment", group = "treatment")) +
        ggplot2::geom_point() +
        ggplot2::geom_line(data = theoretical_missing,
                           ggplot2::aes_string("time", "missing",
                                               color = "treatment", group = "treatment"),
                           linetype = "dashed") +
        ggplot2::labs(title = "Dropout", y = "Proportion dropout", x = "Time point") +
        ggplot2::ylim(0,1) +
       ggplot2::theme_minimal()


}
.plot_trend <- function(paras) {
    time <- get_time_vector(paras)

    ES <- get_effect_size(paras)
    subtitle <- paste("The treatment effect at posttest",
                      ifelse(ES$standardizer == "raw", "", " (Cohen's d)"),
                      " = ",
                      ES$ES,
                      sep = "")

    caption <- ifelse(ES$standardizer == "raw",
                      "N.B.: The treatment effect is the raw (unstandardized) difference",
                      paste("N.B.: Cohen's d is calculated using the", ES$standardizer))

    y <- paras$fixed_intercept + paras$fixed_slope * time
    y1 <-  paras$fixed_intercept + (paras$fixed_slope + get_slope_diff(paras)/paras$T_end) * time

    d <- data.frame(y = c(y, y1),
                    time = rep(time, 2),
                    treatment = rep(c(0, 1), each = length(y)))
    d$treatment <- factor(d$treatment, labels = c("Control", "Treatment"))

    ggplot2::ggplot(d, ggplot2::aes_string("time", "y", color = "treatment")) +
        ggplot2::geom_line(show.legend = TRUE) +
        ggplot2::geom_point(show.legend = FALSE) +
        ggplot2::labs(title = "Treatment effects", y = "Outcome", x = "Time point",
                      subtitle = subtitle,
                      caption = caption) +
        ggplot2::theme_minimal()
}
.plot_link <- function(object, RE_level, show = TRUE, ...) {
    # To get RE intervals
    m <- marginalize(object, link_scale = TRUE, ...)
    p <- plot.plcp_marginal_nested(m,
                                   RE = TRUE,
                                   RE_level = RE_level,
                                   link_scale = TRUE,
                                   ...)
    if(show) {
        plot(p)
    }
    return(invisible(p))
}
.plot_link_ridges <- function(object, RE_level, show = TRUE, ...) {
    # To get RE intervals
    m <- marginalize(object, link_scale = TRUE, ...)
    p <- plot.plcp_marginal_nested(m,
                                   type = "trend_ridges",
                                   RE = TRUE,
                                   RE_level = RE_level,
                                   link_scale = TRUE,
                                   ...)
    if(show) {
        plot(p)
    }
    return(invisible(p))
}

.plot_diff <- function(x, ...) {
    m <- marginalize(x, link_scale = TRUE)
    .plot_diff_marg(m, ...)
}
.plot_diff_marg <- function(x, type = "post_diff", hu = FALSE, fixed_overall = NULL) {

    if(hu) {
        tmp <- x$post_hu_ps$effect
        ES <- x$post_hu[x$post_hu$var == "marg_hu_post_diff", "est"]
        ES_med <- x$post_hu[x$post_hu$var == "median_hu_post_diff", "est"]
        ES_ratio <- x$post_hu[x$post_hu$var == "marg_OR", "est"]
        ES_ratio_med <- x$post_hu[x$post_hu$var == "median_OR", "est"]
    } else {
        tmp <- x$post_ps$effect
        ES <- x$post[x$post$var == "marg_post_diff", "est"]
        ES_med <- x$post[x$post$var == "median_post_diff", "est"]
        if(x$paras$family == "binomial") {
            # Use OR
            ES_ratio <- x$post[x$post$var == "marg_OR", "est"]
            ES_ratio_med <- x$post[x$post$var == "median_OR", "est"]
            tmp$ratio <- tmp$OR
        } else {
            ES_ratio <- x$post[x$post$var == "marg_RR", "est"]
            ES_ratio_med <- x$post[x$post$var == "median_RR", "est"]
        }

    }

    #tmp$fill <- ifelse(tmp$percentile == 0.5, "median", "other")
    #tmp$fill[which.min(abs(tmp$diff - ES))] <- "mean"


    if(type == "post_diff" | type == "post_diff_ratio") {

        if(abs(ES_med - ES) < 0.0001) {
            breaks <- ES
            labels <- paste(round(ES, 2), " \n(mean,\nmedian)")
        } else {
            breaks <- c(ES_med, ES)
            labels <- c(paste(round(ES_med, 2), " (median)"),
                        paste(round(ES, 2), " (mean)")
            )
        }

        p0 <- ggplot(tmp, aes(percentile, diff, fill = fill)) +
            geom_histogram(stat = "identity", color = "white", fill = "#3498db", alpha = .75) +
            #geom_hline(yintercept = 0, linetype = "solid", size = 0.75) +
            geom_hline(yintercept = ES, linetype = "dotted", alpha = 0.75, size = 0.75) +
            geom_hline(yintercept = ES_med, linetype = "dashed", alpha = 0.75, size = 0.75) +
            scale_y_continuous(sec.axis = sec_axis(~ ., breaks = breaks,
                                                   labels = labels
            )
            ) + theme_minimal() +
            theme(legend.position = "none")

        if(!is.null(fixed_overall)) {
            p0 <- p0 + geom_hline(yintercept = fixed_overall$diff, color = "#e74c3c")
        }
    }


    # Ratio
    if(type == "post_ratio" | type == "post_diff_ratio") {
        if(hu) tmp$ratio <- tmp$OR
        if(abs(ES_ratio_med - ES_ratio) < 0.0001) {
            breaks <- ES_ratio
            labels <- paste(round(ES_ratio, 2), " \n(mean,\nmedian)")
        } else {
            breaks <- c(ES_ratio_med, ES_ratio)
            labels <- c(paste(round(ES_ratio_med, 2), " (median)"),
                        paste(round(ES_ratio, 2), " (mean)")
            )
        }

        p1 <- ggplot(tmp, aes(percentile, ratio)) +
            geom_histogram(stat = "identity", color = "white", fill = "#3498db", alpha = .75) +
            #geom_hline(yintercept = 0, linetype = "dotted", size = 0.75) +
            geom_hline(yintercept = ES_ratio, linetype = "dotted", size = 0.75) +
            geom_hline(yintercept = ES_ratio_med, linetype = "dashed", size = 0.75) +
            scale_y_continuous(sec.axis = sec_axis(~ ., breaks = breaks,
                                                   labels = labels
            )) +
            theme_minimal() +
            theme(legend.position = "none")

        if(!is.null(fixed_overall)) {
            p1 <- p1 + geom_hline(yintercept = fixed_overall$ratio, color = "#e74c3c")
        }

    }

    # Return
    if(type == "post_diff") {
        plot(p0)
        return(invisible(list("post_diff" = p0)))
    } else if(type == "post_ratio") {
        plot(p1)
        return(invisible(list("post_ratio" = p1)))
    } else if(type == "post_diff_ratio") {
        gridExtra::grid.arrange(p0, p1)
        return(invisible(list("post_diff" = p0,
                              "post_ratio" = p1)))
    }

}

# Plot design
#' Plot method for \code{study_parameters}-objects
#' @param x An object of class \code{plcp}.
#' @param n specifies which row \code{n} should be used if \code{object}
#' @paran type ...
#' @param fixed_subject_percentiles ...
#' @param fixed_cluster_percentiles ...
#'  is a \code{data.frame} containing multiple setups.
#' @param type indicated what plot to show. If \code{effect} the plot showing the treatment groups
#' change over time will be shown, if \code{dropout} the missing data pattern will be shown,
#' if \code{both} both plots will be shown.
#'
#' @param ... Optional arguments.
#' @export
plot.plcp_nested <- function(x, n = 1, type = "trend", ..., RE = TRUE, RE_level = 2, hu = FALSE) {
    check_installed("ggplot2")
    paras <- x
     if(is.data.frame(paras)) {
          paras <- as.list(paras[n,])
          paras <- do.call(study_parameters, paras)
          #class(paras) <- append(c("plcp"), class(paras))
     }
     if(type == "trend") {
         if(RE) {
             .plot_link(paras,
                        RE_level = RE_level,
                        ...) +
                 labs(y = "y (link scale)")
         } else {
             .plot_trend(paras, ...)
         }
     } else if(type == "trend_ridges") {
         .plot_link_ridges(paras,
                           RE_level = RE_level,
                           ...)

     } else if(type == "dropout") {
         .plot_dropout(paras)
     } else if(type %in% c("post_diff", "post_ratio", "post_ratio_diff")) {
        .plot_diff(x, type = type, hu = hu)
     }

}


.plot_marg <- function(x, Q_long, ymin, ymax, RE = TRUE, overlay = FALSE, ...) {

    x$treatment <- factor(x$treatment, labels = c("Control", "Treatment"))
    Q_long$treatment <- factor(x$treatment, labels = c("Control", "Treatment"))

    if(overlay) {
        # Overlay L1 trajectory on L2 panel
        tmp <- x[x$var == "Within-subject", ]

        # silently ignore overlay when RE_level != 1
        if(nrow(tmp) > 0) {
            tmp$var <- "Subject"
            tmp$color <- "L1"
        }


        x$color <- NA
        x[x$var == "Within-subject", "color"] <- "L1"
        x[x$var == "Subject", "color"] <- "L2"
        x[x$var == "Cluster", "color"] <- "L3"

        x$color <- factor(x$color, levels = c("L1", "L2", "L3"))

        x <- rbind(x, tmp)
    }

    plot_struct <- list(
                        scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")),
                        guides(color = guide_legend(override.aes = list(fill = NA))),
                        scale_fill_brewer(palette = "PuBu"), # PuBu
                        scale_color_manual(values = c("#192a56", "#e84118", "#e84118")),
                        theme_minimal())

    if(RE) {
        ggplot(x, aes(time, mean, group = treatment)) +
            geom_ribbon(data = Q_long, aes(ymin = min,
                                           ymax = max,
                                           y = NULL,
                                           x = time,
                                           group = interaction(width, treatment),
                                           fill = width),
                        alpha = 0.75) +
            geom_line(aes(y = Q50,
                          color = color,
                          linetype = "median",
                          fill = NULL,
                          group = interaction(color, var, treatment)),
                      size = 1) +
            geom_line(aes(color = color,
                          linetype = "mean",
                          fill = NULL,
                          group = interaction(color, var, treatment)),
                      size = 1) +
            # geom_point(aes(y = Q50,
            #                color =
            #                    color)) +
            #scale_color_manual(values = c("median" = "red", "mean" = "red")) +
            facet_wrap(~treatment, ncol = 2) +
            plot_struct
    } else {
        ggplot(x, aes(time, mean, group = treatment, color = treatment)) +
            geom_line(aes(linetype = "mean", fill = NULL), size = 1) +
            geom_line(aes(y = Q50, linetype = "median", fill = NULL), size = 1) +
            geom_point(aes(y = Q50)) +
            plot_struct
    }


}

.make_nested_trend <- function(object, RE, RE_level,
                               var1 = "y",
                               var2 = "y2",
                               var3 = "y3",
                               level1_func = .sample_level1_nested,
                               ...) {
    y1 <- NULL
    y2 <- NULL
    y3 <- NULL

    if(any(RE_level == 1)) {
        x1 <- do.call(level1_func, list(pars = object$paras, ...))
        y1 <- x1[[var1]]
    }
    if(any(RE_level == 2)) {
        y2 <- object[[var2]]
    }
    if(any(RE_level == 3)) {
        y3 <- object[[var3]]

    }
    args <- list("within-subject" = y1,
                 "subject" = y2,
                 "cluster" = y3)
    args <- args[!vapply(args, is.null, logical(1))]
    x <- .rbind_lists(args)
    Q_long <- .rbind_lists(args, func = reshape_eta_sum)

    # Get limits
    if(RE) {
        lims <- .get_facet_lims(d = Q_long,
                                var_names = names(args),
                                min_cols = "min",
                                max_cols = "max")
    } else {
        lims <- .get_facet_lims(d = x,
                                var_names = names(args),
                                min_cols = c("mean", "Q50"),
                                max_cols = c("mean", "Q50"))
    }

    # use same limits for lvl 1, 2 and 3
    if(all(c("subject", "cluster") %in% lims$var)) {
        tmp <- lims[lims$var %in% c("within-subject", "subject", "cluster"), ]
        lims[lims$var == "within-subject", "mean"] <- c(min(tmp$mean), max(tmp$mean))
        lims[lims$var == "subject", "mean"] <- c(min(tmp$mean), max(tmp$mean))
        lims[lims$var == "cluster", "mean"] <- c(min(tmp$mean), max(tmp$mean))
    }

    x$var <- factor(x$var,
                    labels = c("Within-subject", "Subject", "Cluster"),
                    levels = c("within-subject","subject", "cluster"))
    Q_long$var <- factor(Q_long$var,
                         labels = c("Within-subject", "Subject", "Cluster"),
                         levels =  c("within-subject","subject", "cluster"))

    x$color <- x$var
    Q_long$color <- Q_long$var
    lims2 <- lims
    lims2$treatment <- "Control"
    lims <- rbind(lims, lims2)

    lims$var <- factor(lims$var, labels = c("Within-subject", "Subject", "Cluster"),
                       levels = c("within-subject","subject", "cluster"))
    lims$treatment <- factor(lims$treatment, labels = c("Control", "Treatment"))
    x$treatment <- factor(x$treatment, labels = c("Control", "Treatment"))
    Q_long$treatment <- factor(Q_long$treatment, labels = c("Control", "Treatment"))

    list(x = x,
         Q_long = Q_long,
         lims = lims)
}

.plot_nested_trend_ridges <- function(res, trend, RE, RE_level, stat = "density", family, link_scale = FALSE, ...) {

    lims <- trend$lims
    tmp <- lims[lims$var %in% c("Within-subject", "Subject", "Cluster"), ]
    lims[lims$var == "Within-subject", "mean"] <- c(min(tmp$mean), max(tmp$mean))
    lims[lims$var == "Subject", "mean"] <- c(min(tmp$mean), max(tmp$mean))
    lims[lims$var == "Cluster", "mean"] <- c(min(tmp$mean), max(tmp$mean))

    if(family != "binomial") res <- subset(res, y > lims$mean[1] & y < lims$mean[2])

    p <- ggplot(res, aes(x = y, y = time, group = interaction(time, treatment, var), fill = treatment, color = treatment))

    if(family == "binomial" & !link_scale) {
        p <- p +
            ggridges::geom_density_ridges(data = subset(res, var == "Within-subject"),
                                          scale = 0.7, stat = "binline",
                                          aes(height = ..count..),
                                          bins = 20,
                                          rel_min_height = 0.01,
                                          color = alpha("white", 0.33),
                                          alpha = 0.75,
                                          size = 0.3)

    } else if(family == "poisson") {
        p <- p +
            ggridges::geom_density_ridges(data = subset(res, var == "Within-subject"),
                                          scale = 0.7, stat = "binline",
                                          aes(height = ..count..),
                                          binwidth = 1,
                                          rel_min_height = 0.01,
                                          color = alpha("white", 0.33),
                                          alpha = 0.75,
                                          size = 0.3)

    }
    else {
        p <- p +
            ggridges::geom_density_ridges(data = subset(res, var == "Within-subject"),
                                          scale = 0.7, stat = stat,
                                          aes(height = ..count..),
                                          rel_min_height = 0.01,
                                          color = alpha("white", 0.33),
                                          alpha = 0.75,
                                          size = 0.3)
    }

    p + ggridges::geom_density_ridges(data = subset(res, var  %in% c("Subject", "Cluster")), scale = 0.7, stat = stat,
                                      aes(height = ..count..),
                                      #binwidth = 1,
                                      rel_min_height = 0.01,
                                      color = alpha("white", 0.33),
                                      alpha = 0.75,
                                      size = 0.3) +
        geom_path(data = trend$x,
                  aes(x = mean,
                      y = time,
                      linetype = "mean",
                      fill = NULL,
                      group = interaction(treatment, var)),
                  size = 1) +
        geom_path(data = subset(trend$x, var %in% c("Subject", "Cluster")),
                  aes(x = Q50,
                      y = time,
                      linetype = "median",
                      fill = NULL,
                      group = interaction(treatment, var)), size = 1) +
        geom_blank(data = lims, aes(x = mean, y=time)) +
        coord_flip() +
        theme_minimal() +
        #scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(breaks = unique(sort(res$time))) +
        facet_wrap(~var, ncol = 2) +
        scale_fill_manual(values = c("#30394F", "#6ACEEB")) +
        scale_color_manual(values = c("#30394F", "#c0392b"))

}

plot.plcp_marginal_nested <- function(object, type = "trend", ..., RE = TRUE, RE_level = 2, hu = FALSE) {
    check_installed("ggplot2")

    ## DROPOUT
    if(type == "dropout") {
        .plot_dropout(object$paras)
    ## TREND
    } else if(type == "trend") {
        trend <- .make_nested_trend(object = object,
                                    RE = RE,
                                    RE_level = RE_level,
                                    ...)

        if(RE) {
            facets <- list(facet_wrap(var ~ treatment, scales = "free", ncol = 2))
        } else {
            facets <- list(facet_wrap(~var, ncol = 1, scales = "free"))
        }
        p <- .plot_marg(x = trend$x,
                        Q_long = trend$Q_long,
                        RE = RE,
                        ymin = NA,
                        ymax = NA,
                        ...) +
            geom_blank(data = trend$lims) +
            labs(linetype = "",
                 color = "",
                 y = "Y",
                 title = "Change over time") +
            facets

        plot(p)
        return(invisible(p))
    ## POST
    } else if(type %in% c("post_diff", "post_ratio", "post_ratio_diff")) {
        .plot_diff_marg(object, type = type, ...)

        ## Ridges
    } else if(type == "trend_ridges") {
        check_installed("ggridges")
        res <- .mu_vec_to_long(object,
                               RE_level = RE_level,
                               ...)
        trend <- .make_nested_trend(object = object,
                                    RE = RE,
                                    RE_level = RE_level,
                                    ...)
        .plot_nested_trend_ridges(res = res,
                                  trend = trend,
                                  RE_level = RE_level,
                                  family = object$paras$family,
                                  ...)
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
