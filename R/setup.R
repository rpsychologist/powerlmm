#' Setup study parameters
#'
#' Setup the parameters for calculating power for longitudinal multilevel studies
#' comparing two groups. Ordinary two-level models (subjects with repeated measures),
#' and longitudinal three-level models with clustering due to therapists, schools, provider etc,
#' are supported. Random slopes at the subject level and cluster level are
#' possible. Cluster sizes can be unbalanced, and vary by treatment.
#' Partially-nested designs are supported. Missing data can also be accounted
#' for.
#'
#' @param n1 Number of level 1 units, e.g. measurements per subject.
#' @param n2 Number of level 2 units per level 3 unit, e.g. subjects per cluster.
#' Unbalanced cluster sizes are supported, see \code{\link{unequal_clusters}}.
#' @param n3 Number of level 3 units per treatment, can be different in each
#'  treatment arm, see \code{\link{per_treatment}}.
#' @param T_end Time point of the last measurement. If \code{NULL} it will be set
#' to \code{n1 - 1}.
#' @param fixed_intercept Average baseline value, assumed to be equal for both groups.
#' @param fixed_slope Overall change per unit time, in the control group.
#' @param sigma_subject_intercept Subject-level random intercept.
#' @param sigma_subject_slope Subject-level random slope.
#' @param sigma_cluster_intercept Cluster-level random intercept.
#' @param sigma_cluster_slope Cluster-level random slope.
#' @param sigma_error Within-subjects (residual) variation.
#' @param icc_slope Proportion of slope variance
#' at the cluster level.
#' @param var_ratio Ratio of the random
#' slope variance to the within-subject variance.
#' @param icc_pre_subject Amount of baseline
#' variance at the subject level. N.B. this is not really an ICC.
#' @param icc_pre_cluster Amount of baseline
#' variance at the cluster level.
#' @param cor_subject Correlation between the subject-level random intercept
#'  and slopes.
#' @param cor_cluster Correlation between the cluster-level random intercept
#' and slopes.
#' @param cor_within Correlation of the level 1 residual. Currently ignored in
#' the analytical power calculations.
#' @param cohend Standardized between-group treatment effect at the last time point,
#' see \emph{Details}.
#' @param partially_nested \code{logical}; indicates if there's clustering in both
#' arms or only in the treatment arm.
#' @param dropout Dropout process, see \code{\link{dropout_weibull}} or
#' \code{\link{dropout_manual}}. Assumed to be 0 if \code{NULL}.
#' @param deterministic_dropout \code{logical}; if \code{FALSE} the input to
#' \code{dropout} will be treated as stochastic and dropout will sampled
#' from a multinomial distribution.
#' @return A \code{list} or \code{data.frame} of parameters values, either of
#' class \code{plcp} or \code{plcp_multi} if multiple parameters are compared.
#'
#' @details
#'
#' \bold{Comparing a combination of parameter values}
#'
#' It is possible to setup a grid of parameter combinations by entering the values
#' as vectors. All unique combinations of the inputs will be returned. This is
#' useful if you want see how different values of the parameters affect power.
#' See also the convenience function \code{\link{get_power_table}}.
#'
#' \bold{Standardized and unstandardized inputs}
#'
#' All parameters of the models can be specified. However, many of the raw
#' parameter values in a multilevel/LMM do no directly affect the power of the
#' test of the \code{treatment:time}-coefficient. Power will depend greatly on the relative
#' size of the parameters, therefore, it is possible to setup your calculations
#' using only standardized inputs, or by a combination of raw inputs and
#' standardized inputs. For instance, if \code{sigma_subject_slope} and
#' \code{icc_slope} is specified, the \code{sigma_cluster_slope} will be
#' solved for. More information regarding the standardized inputs are available in
#' the two-level and three-level vignettes.
#'
#' \bold{Cohen's d calculation}
#'
#' Cohen's \emph{d} is calculated by using the baseline standard deviation as the denominator.
#' The choice of denominator differs between fields, and other options will be added in
#' future releases.
#'
#' \bold{Two- or three-level models}
#'
#' If either \code{sigma_cluster_slope} or \code{icc_slope} and
#'  \code{sigma_cluster_intercept} or \code{icc_pre_cluster} is zero or
#' \code{NULL} it will be assumed a two-level design is wanted.
#'
#' \bold{Unequal cluster sizes and unbalanced allocation}
#'
#' It is possible to specify different cluster sizes using
#' \code{\link{unequal_clusters}}. Cluster sizes can vary between treatment arms
#' by also using \code{\link{per_treatment}}. The number of clusters per treatment can
#' also be set by using \code{\link{per_treatment}}. See \code{\link{per_treatment}}
#' and \code{\link{unequal_clusters}} examples of their use.
#'
#' \bold{Missing data and dropout}
#'
#' Accounting for missing data in the power calculations is possible. Currently,
#' \code{dropout} can be specified using either \code{\link{dropout_weibull}} or
#'  \code{\link{dropout_manual}}. It is possible to have different dropout
#'  patterns per treatment group using \code{\link{per_treatment}}. See their
#'  respective help pages for examples of their use.
#'
#' @seealso \code{\link{get_power}}, \code{\link{simulate.plcp}}
#'
#'
#' @examples
#' # Three level model with both subject- and cluster-level random slope
#' # Power calculation using standardized inputs
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 4,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)
#'
#' get_power(p)
#'
#' # The same calculation with all parameters specified directly
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 4,
#'                       T_end = 10,
#'                       fixed_intercept = 37,
#'                       fixed_slope = -0.65,
#'                       sigma_subject_intercept = 2.8,
#'                       sigma_subject_slope = 0.4726944,
#'                       sigma_cluster_intercept = 0,
#'                       sigma_cluster_slope = 0.1084435,
#'                       sigma_error = 2.8,
#'                       cor_subject = -0.5,
#'                       cor_cluster = 0,
#'                       cohend = -0.8)
#' get_power(p)
#'
#' # Standardized and unstandardized inputs
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 4,
#'                       sigma_subject_intercept = 2.8,
#'                       icc_pre_cluster = 0.07,
#'                       sigma_subject_slope = 0.47,
#'                       icc_slope = 0.05,
#'                       sigma_error = 2.8,
#'                       cohend = -0.8)
#'
#' get_power(p)
#'
#' ## Two-level model with subject-level random slope
#' p <- study_parameters(n1 = 11,
#'                       n2 = 40,
#'                       icc_pre_subject = 0.5,
#'                       var_ratio = 0.03,
#'                       cohend = -0.8)
#' get_power(p)
#'
#' # add missing data
#' p <- update(p, dropout = dropout_weibull(0.2, 1))
#' get_power(p)
#'
#' ## Comparing a combination of values
#' p <- study_parameters(n1 = 11,
#'                       n2 = c(5, 10),
#'                       n3 = c(2, 4),
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = c(0, 0.05),
#'                       cohend = c(-0.5, -0.8))
#'
#' get_power(p)
#' @export
study_parameters <- function(n1, n2, n3=1, T_end=NULL,
                             fixed_intercept = 0L,
                             fixed_slope = 0L,
                             sigma_subject_intercept = NULL,
                             sigma_subject_slope = NULL,
                             sigma_cluster_intercept = NULL,
                             sigma_cluster_slope = NULL,
                             sigma_error = NULL,
                             cor_subject = 0L,
                             cor_cluster = 0L,
                             cor_within = 0L,
                             var_ratio = NULL,
                             icc_slope = NULL,
                             icc_pre_subject = NULL,
                             icc_pre_cluster = NULL,
                             cohend = 0L,
                             partially_nested = FALSE,
                             dropout = 0L,
                             deterministic_dropout = TRUE) {

    #if(!is.per_treatment(n2) & length(n2) == 1) n2 <- list(n2)

    # drop out checks
    if(is.numeric(dropout) && any(dropout != 0)) stop("'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'")
    if(is.per_treatment(dropout)) {
        tx <- dropout[[1]]$treatment
        cc <- dropout[[1]]$control
        if(is.numeric(cc) && any(cc != 0)) stop("Control group's 'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'")
        if(is.numeric(tx) && any(tx != 0)) stop("Treatment group's 'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'")
    }

    # warn n3 is ignores
    if(is.unequal_clusters(n2) & !is.per_treatment(n2) & is.per_treatment(n3)) {
        message("'n3' per_treatment argument is ignored. 'n3' is automatically based on length of unequal_clusters")
    }

    # Check if enough parameters are provided
    if(!is.null(icc_slope) & (!is.null(sigma_subject_slope) &
                              !is.null(sigma_cluster_slope))) {
        stop("Can't use 'icc_slope' with both 'sigma_subject_slope' and 'sigma_cluster_slope'")

    }
    if(!is.null(sigma_subject_slope) & !is.null(icc_pre_subject) & !is.null(var_ratio)) {
        if(any(var_ratio == 0)) stop("'var_ratio' can't be zero")
    }
    if(!is.null(sigma_subject_slope) &
       !is.null(var_ratio) &
       !is.null(sigma_error)) {
        stop("'sigma_subject_slope' or 'var_ratio' or 'sigma_error' should be NULL")

    }

    if(!is.null(sigma_subject_intercept) & !is.null(icc_pre_subject)) {
        stop("Can't use both 'icc_pre_subject' and 'sigma_subject_intercept'")
    }

    if(!is.null(sigma_cluster_intercept) & !is.null(icc_pre_cluster)) {
        stop("Can't use both 'icc_pre_cluster' and 'sigma_cluster_intercept'")
    }

    if(is.null(sigma_subject_intercept) & is.null(icc_pre_subject)) {
         stop("Both 'sigma_subject_intercept' and 'icc_pre_subject' can't be NULL")
    }
    if(is.null(sigma_subject_slope) & (is.null(var_ratio) || var_ratio == 0) &
       is.null(sigma_cluster_slope) & is.null(icc_slope) &
       is.null(icc_pre_subject) & is.null(icc_pre_cluster)) {
        stop("Both 'sigma_subject_slope' and 'var_ratio' can't be NULL")
    }
    if(!(is.null(icc_slope) || icc_slope == 0) & (is.null(var_ratio) || var_ratio == 0) &
       (is.null(sigma_cluster_slope) & is.null(sigma_subject_slope))) {
        stop("Argument 'icc_slope' requires that 'var_ratio' is specified and > 0.")
    }
    if(!is.null(icc_pre_cluster) & !is.null(icc_pre_subject)) {
        if(any((icc_pre_subject + icc_pre_cluster) >= 1)) {
            stop("'icc_pre_subject' and 'icc_pre_cluster' shouldn't sum to >= 1")
        }
    }
    # if(!is.null(icc_pre_subject) & !is.null(icc_pre_cluster) & !is.null(var_ratio) &
    #    is.null(icc_slope) & is.null(sigma_subject_slope) & is.null(sigma_cluster_slope) &
    #    is.null(sigma_error)) {
    #     stop("Not enough information provided.")
    # }

    no_raw <- is.null(sigma_subject_intercept) &
                is.null(sigma_subject_slope) &
                is.null(sigma_cluster_intercept) &
                is.null(sigma_cluster_slope) &
                is.null(sigma_error)

    if(!(is.null(icc_slope) || icc_slope == 0) &
       !(is.null(var_ratio) || var_ratio == 0) &
       !(no_raw) &
       (is.null(icc_pre_subject) | is.null(icc_pre_cluster)) &
       is.null(sigma_error) &
       is.null(sigma_subject_slope) &
       is.null(sigma_cluster_slope)) {
        stop("Either 'sigma_error', 'sigma_subject_slope', or 'sigma_cluster_slope' need to be specified")
    }
    if(!is.null(icc_slope) & (!is.null(sigma_subject_slope) | !is.null(sigma_cluster_slope)) &
       (is.null(sigma_error) || sigma_error == 0) & (is.null(var_ratio) || var_ratio == 0)) {
        stop("'sigma_error' or 'var_ratio' should be specified and > 0")
    }


    #if(is.null(icc_pre_cluster) || all(icc_pre_cluster == 0)) icc_pre_cluster <- NULL
    #if(is.null(icc_pre_subject) || all(icc_pre_subject == 0)) icc_pre_subject <- NULL

    args <- list(
        n1 = n1,
        n2 = n2,
        n3 = n3,
        T_end = T_end,
        fixed_intercept = fixed_intercept,
        fixed_slope = fixed_slope,
        sigma_subject_intercept = sigma_subject_intercept,
        sigma_subject_slope = sigma_subject_slope,
        sigma_cluster_intercept = sigma_cluster_intercept,
        sigma_cluster_slope = sigma_cluster_slope,
        sigma_error = sigma_error,
        cor_subject = cor_subject,
        cor_cluster = cor_cluster,
        cor_within = cor_within,
        var_ratio = var_ratio,
        icc_pre_subject = icc_pre_subject,
        icc_pre_cluster = icc_pre_cluster,
        icc_slope = icc_slope,
        cohend = cohend,
        partially_nested = partially_nested,
        dropout = dropout,
        deterministic_dropout = deterministic_dropout
    )
    save_call <- args

    #args$dropout <- lapply(args$dropout, make_list_weibull)
    # if(is.list(args$dropout) & !is.per_treatment(args$dropout)) {
    #     args$dropout <- list(dropout)
    # }

#    unequal_clust <- lapply(seq_along(n2), function(i) is.unequal_clusters(n2[i]))
 #   unequal_clust <- unlist(unequal_clust)

  #  if(any(unequal_clust) & length(n3) > 1) stop("Can't combine `unequal_clusters` with different `n3` values.")
    tmp_args <- args[!vapply(args, is.null, logical(1))]

    tmp <- expand.grid(tmp_args)
    if(is.null(args$T_end)) tmp$T_end <- tmp$n1 - 1

    # two level
    if(!is.null(sigma_subject_slope) &
       !(is.null(var_ratio) || var_ratio == 0) &
       is.null(sigma_error) &
       (is.null(icc_slope) || icc_slope == 0)) {
        x <- with(tmp, sqrt(sigma_subject_slope^2 / var_ratio))
        sigma_error <- x
        tmp$sigma_error <- x
    }
    if(!is.null(icc_pre_subject) &
       (is.null(icc_pre_cluster) || icc_pre_cluster == 0) &
       is.null(sigma_error) &
       (is.null(icc_slope) || icc_slope == 0)) {

        tmp$sigma_error <- sqrt(1 - (tmp$icc_pre_subject))
        sigma_error <- tmp$sigma_error
        tmp$sigma_subject_intercept <- sqrt(tmp$icc_pre_subject)
        tmp$sigma_cluster_intercept <- 0L
    }


    # two level

    if(is.null(sigma_cluster_intercept) & (is.null(icc_pre_cluster) || icc_pre_cluster == 0)) {
        tmp$sigma_cluster_intercept <- 0L
    }
    if(is.null(sigma_cluster_slope) & (is.null(icc_slope) || icc_slope == 0)) {
        tmp$sigma_cluster_slope <- 0L
    }
    if(!is.null(icc_pre_subject) &
       !is.null(sigma_error) &
       is.null(sigma_subject_intercept) &
       is.null(sigma_cluster_intercept) &
       (is.null(icc_pre_cluster) || icc_pre_cluster == 0)) {

        tmp$sigma_subject_intercept <- with(tmp,
                                            (icc_pre_subject * sigma_error^2) /
                                                (1-icc_pre_subject))
        tmp$sigma_subject_intercept <- sqrt(tmp$sigma_subject_intercept)

    }
    if(!is.null(icc_pre_subject) &
       !is.null(var_ratio) &
       !is.null(sigma_subject_slope) &
       is.null(sigma_error) &
       is.null(sigma_subject_intercept) &
       is.null(sigma_cluster_intercept) &
       (is.null(icc_pre_cluster) || icc_pre_cluster == 0)) {

        tmp$sigma_subject_intercept <- with(tmp,
                                            (icc_pre_subject * sigma_error^2) /
                                                (1-icc_pre_subject))
        tmp$sigma_subject_intercept <- sqrt(tmp$sigma_subject_intercept)

    }
    if(!is.null(var_ratio) &
       (is.null(icc_slope) || icc_slope == 0) &
       !is.null(sigma_error)) {
        tmp$sigma_subject_slope <- sqrt(with(tmp, var_ratio *
                                                 sigma_error^2))
    }



    # three level

    if(!is.null(icc_pre_subject) &
       !is.null(icc_pre_cluster) &
       is.null(sigma_error)) {
        tmp$sigma_error <- sqrt(1-(tmp$icc_pre_subject +
                                       tmp$icc_pre_cluster))
        tmp$sigma_subject_intercept <- sqrt(tmp$icc_pre_subject)
        tmp$sigma_cluster_intercept <- sqrt(tmp$icc_pre_cluster)


    }
    if(!is.null(var_ratio) & !is.null(icc_slope)) {
        if(is.null(sigma_error) & !is.null(sigma_cluster_slope)) {
            sigma_error <- sigma_cluster_slope^2/(icc_slope * var_ratio)
            sigma_error <- sqrt(sigma_error)
            tmp$sigma_error <- sigma_error
        }
        if(is.null(sigma_error) & !is.null(sigma_subject_slope)) {
            sigma_error <- sigma_subject_slope^2/((1-icc_slope) * var_ratio)
            sigma_error <- sqrt(sigma_error)
            tmp$sigma_error <- sigma_error
        }

        v1 <- with(tmp, var_ratio * sigma_error^2 * icc_slope)
        tmp$sigma_cluster_slope <- sqrt(v1)
        tmp$sigma_subject_slope <- sqrt(with(tmp, var_ratio *
                                                 sigma_error^2 * (1-icc_slope)))

    }
    if(is.null(var_ratio) & !is.null(icc_slope) & !is.null(sigma_error)) {
        if(!is.null(sigma_cluster_slope)) {
            x <- with(tmp, sigma_cluster_slope^2/icc_slope)
            u1 <- x - tmp$sigma_cluster_slope^2
            u1 <- sqrt(u1)
            tmp$sigma_subject_slope <- u1
        }
        if(!is.null(sigma_subject_slope)) {
            x <- with(tmp, sigma_subject_slope^2/(1-icc_slope))
            v1 <- x - tmp$sigma_subject_slope^2
            v1 <- sqrt(v1)
            tmp$sigma_cluster_slope <- v1
        }

    }

    # no lvl 3 random slope
    if(!is.null(var_ratio) &
       (is.null(icc_slope) || icc_slope == 0) &
       is.null(sigma_error) &
       is.null(sigma_cluster_slope) & is.null(sigma_subject_slope)) {
        tmp$sigma_subject_slope <- sqrt(with(tmp, var_ratio *
                                                 sigma_error^2))


    }
    # var ratio also null
    if((all(is.null(var_ratio)) || all(var_ratio == 0)) &
      (is.null(icc_slope) || icc_slope == 0) &
       is.null(sigma_cluster_slope) & is.null(sigma_subject_slope)) {
        tmp$sigma_subject_slope <- 0L
    }


    if(!is.null(icc_pre_subject) & !is.null(icc_pre_cluster)) {
        tmp$sigma_subject_intercept <-
            sqrt(
                tmp$sigma_error^2 /
                    (1-(tmp$icc_pre_subject + tmp$icc_pre_cluster)) *
                    tmp$icc_pre_subject
            )
        tmp$sigma_cluster_intercept <-
            sqrt(
                tmp$sigma_error^2 /
                    (1-(tmp$icc_pre_subject + tmp$icc_pre_cluster)) *
                    tmp$icc_pre_cluster
            )
    }
    if(is.null(icc_pre_subject) & !is.null(icc_pre_cluster) &
       !is.null(sigma_subject_intercept)) {
           v0 <- with(tmp, (icc_pre_cluster*(sigma_error^2 + sigma_subject_intercept^2))
                 /(1-icc_pre_cluster))
           tmp$sigma_cluster_intercept <- sqrt(v0)
    }
    if(is.null(icc_pre_cluster) & !is.null(icc_pre_subject) &
       !is.null(sigma_cluster_intercept)) {
        u0 <- with(tmp, (icc_pre_subject*(sigma_error^2 + sigma_cluster_intercept^2))
                   /(1-icc_pre_subject))
        tmp$sigma_subject_intercept <- sqrt(u0)
    }

    cols <- which(colnames(tmp) %in% c("icc_slope",
                                       "var_ratio",
                                       "icc_pre_cluster",
                                       "icc_pre_subject"))
    cols <- colnames(tmp)[cols]
    paras <- tmp[, !(names(tmp) %in% cols)]

    if((is.data.frame(paras) & nrow(paras) == 1)) {
        paras <- as.list(paras)
    }
    if(is.data.frame(paras)) {
        class(paras) <- append(c("plcp_multi"), class(paras))
    } else class(paras) <- append(c("plcp"), class(paras))
    if(is.null(paras$cor_cluster)) paras$cor_cluster <- cor_cluster
    if(is.null(paras$cor_subject)) paras$cor_subject <- cor_subject
    if(all(paras$sigma_cluster_slope == 0) &
       all(paras$sigma_cluster_intercept == 0)) {
        class(paras) <- append(class(paras), c("plcp_2lvl"))
    } else {
        class(paras) <- append(class(paras), c("plcp_3lvl"))
    }

    attr(paras, "call") <- save_call
    paras

}

sim_parameters <- function(...) {
    dots <- list(...)
    warning("sim_parameters is deprecated")
    do.call(study_parameters, dots)
}


print_per_treatment <- function(n, width = 0) {
    x <- lapply(seq_along(n), print_per_treatment_, x = n)
    x <- format(x, width = width)
    x <- paste(x, " (", names(n), ")", sep ="")
    x <- paste(unlist(x), collapse = "\n                   ")
    x
}
print_per_treatment_ <- function(i, x) {
    name <- names(x)[i]
    x <- x[[i]]
    paste(paste(unlist(x), collapse = ", "), sep ="")
}

deparse_n2 <- function(n2) {
    n2_attr <- attr(n2, "func")
    if(!is.null(n2_attr) && (n2_attr != "manual")) {
        n2 <- deparse(n2_attr)
    }
    n2
}
prepare_print_plcp <- function(x, two_level = FALSE) {
    n1 <- x$n1
    n2 <- get_n2(x)
    n2$treatment <- deparse_n2(n2$treatment)
    n2$control <- deparse_n2(n2$control)
    n3 <- get_n3(x)
    tot_n <- get_tot_n(x)
    width <- max(nchar(print_per_treatment_(1, n2)),
                 nchar(print_per_treatment_(2, n2)),
                 nchar(print_per_treatment_(3, tot_n)))
    if(two_level) width <- max(vapply(tot_n, nchar, numeric(1)))
    n2 <- print_per_treatment(n2, width = width)
    n3 <- print_per_treatment(n3, width = width)

    tot_n <- print_per_treatment(tot_n, width = width)

    icc_slope <- round(get_ICC_slope(x), 2)
    var_ratio <- round(get_var_ratio(x), 2)
    icc_pre_clusters <- round(get_ICC_pre_clusters(x), 2)
    icc_pre_subjects <- round(get_ICC_pre_subjects(x), 2)

    cohend <- x$cohend

    gd <- get_dropout(x)
    gd$time <- format(gd$time, nsmall = 0, digits = 3, width = 2)
    gd$control <-  format(gd$control*100, nsmall = 0, digits = 0, width = 2)
    gd$treatment <- format(gd$treatment*100, nsmall = 0, digits = 0, width = 2)
    colnames(gd) <- c("time", "%, control", "%, treatment")
    gd <- print_per_treatment(gd)

    res <- structure(list(n1 = n1,
                          n2 = n2,
                          n3 = n3,
                          total_n = tot_n,
                          dropout = gd,
                          icc_pre_subjects = icc_pre_subjects,
                          icc_pre_clusters = icc_pre_clusters,
                          icc_slope = icc_slope,
                          var_ratio = var_ratio,
                          cohend = cohend,
                          method = "Study setup (three-level)"),
                     class = "power.htest")
    res
}
prepare_print_plcp_2lvl <- function(x) {
    res <- prepare_print_plcp(x, two_level = TRUE)
    if(!is.list(x$dropout)) res$dropout <- "No missing data"
    res$method <- "Study setup (two-level)"
    res$icc_slope <- NULL
    res$icc_pre_clusters <- NULL
    res$n2 <- res$total_n
    res$n3 <- NULL
    res$total_n <- NULL

    res
}
prepare_print_plcp_3lvl <- function(x) {
    res <- prepare_print_plcp(x)
    if(!is.list(x$dropout)) res$dropout <- "No missing data"
    if(x$partially_nested) res$method <- "Study setup (three-level, partially nested)"
    if(is.unequal_clusters(x$n2)) {
        if(is.per_treatment(x$n2)) {
            if(class(x$n2[[1]]$treatment[[1]]) == "plcp_unequal_clusters") {
                n2 <- x$n2[[1]]$treatment[[1]]()
                n2_attr <- attr(n2, "func")
            } else {
                n2 <- x$n2[[1]]$control[[1]]()
                n2_attr <- attr(n2, "func")
            }
        } else {
            n2 <- x$n2[[1]]()
            n2_attr <- attr(n2, "func")
        }

        if(!is.null(n2_attr) & n2_attr != "manual") res$note <- "n2 is randomly sampled"

    }


    res

}

#' Print method for three-level \code{study_parameters}-objects
#' @param x An object of class \code{plcp_3lvl}.
#' @param ... Optional arguments.
#' @method print plcp_3lvl
#' @export
print.plcp_3lvl <- function(x, ...) {
    res <- prepare_print_plcp_3lvl(x)
    print(res, digits = 2)
}

#' Print method for two-level \code{study_parameters}-objects
#' @param x An object of class \code{plcp_2lvl}.
#' @param ... Optional arguments.
#' @method print plcp_2lvl
#' @export
print.plcp_2lvl <- function(x, ...) {
   res <- prepare_print_plcp_2lvl(x)

    print(res, digits = 2, ...)
}

get_slope_diff <- function(paras) {
    with(paras, cohend * sqrt(sigma_subject_intercept^2 + sigma_cluster_intercept^2 + sigma_error^2))
}


# print multi-sim ---------------------------------------------------------

replace_repeating <- function(x, empty) {
    lagx <- x[seq_len(length(x) - 1)]
    lagx <- c(NA, lagx)
    x[x == lagx] <- empty

    x
}
get_dropout_post <- function(object) {
    x <- get_dropout(object)
    x[nrow(x),]
}
prepare_multi_setup <- function(object, empty = ".", digits = 2) {
    paras <- object

    n2 <- lapply(1:nrow(paras), function(i) {
        x <- get_n2(as.plcp(paras[i,]))
        data.frame(treatment = paste(x$treatment, collapse = ","), control = paste(x$control, collapse = ","), stringsAsFactors = FALSE)
    })
    n2 <- do.call(rbind, n2)

    n3 <- lapply(1:nrow(paras), function(i) {
        x <- get_n3(as.plcp(paras[i,]))
    })
    n3 <- do.call(rbind, n3)



    object$icc_pre_cluster <- get_ICC_pre_clusters(object)
    object$icc_pre_subject <- get_ICC_pre_subjects(object)
    object$icc_slope <- get_ICC_slope(object)
    object$var_ratio <- get_var_ratio(object)

    out <- object
    out$icc_pre_cluster <- round(object$icc_pre_cluster, digits)
    out$icc_pre_subject <- round(object$icc_pre_subject, digits)
    out$icc_slope <- round(object$icc_slope, digits)
    out$var_ratio <- round(object$var_ratio, digits)

    dropout <- lapply(1:nrow(object), function(i) get_dropout_post(object[i, ]))
    dropout <- do.call(rbind, dropout)
    if(all(dropout$control == dropout$treatment)) {
        out$dropout <- dropout$treatment
    } else {
        out$dropout_tx <- dropout$treatment
        out$dropout_cc <- dropout$control
        out <- subset(out, select = -dropout)
    }
    if(all(n2$treatment == n2$control)) {
        out$n2 <- n2$treatment
    } else {
        out$n2_tx <- n2$treatment
        out$n2_cc <- n2$control
        out <- subset(out, select = -n2)
    }

    unequal_clust <- lapply(seq_along(object$n2), function(i) is.unequal_clusters(object$n2[i]))
    unequal_clust <- unlist(unequal_clust)

    if(any(unequal_clust)) {
        out <- subset(out, select = -n3)
    } else {
        if(all(n3$treatment == n3$control)) {
            out$n3 <- n3$treatment
        } else {
            out$n3_tx <- n3$treatment
            out$n3_cc <- n3$control
            out <- subset(out, select = -n3)
        }
    }

    out_dense <- out
    for(i in 1:ncol(out)) {
        out[,i] <- replace_repeating(out[,i], empty = empty)
    }

    list(out = out,
         out_dense = out_dense,
         object = object)
}

# prepare
get_multi_title <- function(object) {
    UseMethod("get_multi_title")
}
get_multi_title.plcp_2lvl <- function(object) {
    "# Multi-study setup (two-level)"
}
get_multi_title.plcp_3lvl <- function(object) {
   "# Multi-study setup (three-level)"
}

select_setup_cols <- function(x) {
    cols <- c("n1",
              "n2", "n2_tx", "n2_cc",
              "n3", "n3_tx", "n3_cc",
              "dropout", "dropout_tx", "dropout_cc",
              "icc_pre_subject", "icc_pre_cluster", "icc_slope", "var_ratio", "cohend")
    cols[cols %in% colnames(x)]
}

#' Print method for \code{study_parameters}-multiobjects
#' @param x An object of class \code{plcp_multi}.
#' @param print_max The number of rows to show
#' @param empty Symbol used to replace repeating non-unique parameters
#' @param digits Digits to show
#' @param ... Optional arguments.
#' @method print plcp_multi
#' @export
print.plcp_multi <- function(x, print_max = 10, empty = ".", digits = 2, ...) {
    nr <- nrow(x)
    if(nr <= print_max) rmax <- nr else rmax <- print_max
    hidden_row <- nr - print_max
    x <- x[1:rmax, ]
    pp <- prepare_multi_setup(x, empty = empty, digits = digits)
    out <- pp$out
    out <- as.data.frame(out)
    cat(get_multi_title(pp$object), "\n")

    print(out[, select_setup_cols(out)])
    if(hidden_row > 0) {
        cat("# ...", hidden_row, "setups not shown.")
    }

    invisible(x)
}




# helpers -----------------------------------------------------------------

eval_n2 <- function(n2) {
    n2 <- n2[[1]]()
    n2 <- round(n2, 0)
    func <- attr(n2, "func")
    if(func != "manual") {
        trunc <- attr(n2, "trunc")
        repl <- attr(n2, "replace")

        n2[n2 < trunc] <- repl
    }
    n2 <- n2[n2 > 0]
    if(length(n2) < 1) stop("All clusters of size 0")
    attr(n2, "func") <- func
    n2
}

prepare_paras <- function(paras) {
    paras_tx <- paras
    if (is.per_treatment(paras$n3)) {
        n3_tx <- paras$n3[[1]]$treatment
        n3_cc <- paras$n3[[1]]$control
        paras$n3 <- n3_cc
        paras_tx$n3 <- n3_tx
    }

    if(is.per_treatment(paras$n2)) {

        paras_tx$n2 <- paras$n2[[1]]$treatment
        paras$n2 <- paras$n2[[1]]$control


        if(is.unequal_clusters(paras$n2)) {
            paras$n2 <- eval_n2(paras$n2)
            paras$n3 <- length(paras$n2)

        }
        if(is.unequal_clusters(paras_tx$n2)) {
            paras_tx$n2 <- eval_n2(paras_tx$n2)
            paras_tx$n3 <- length(paras_tx$n2)

        }
    }
    if(is.unequal_clusters(paras$n2)) {
        paras$n2 <- eval_n2(paras$n2)
        paras_tx$n2 <- paras$n2
        paras$n3 <- length(unlist(paras$n2))
        paras_tx$n3 <-  paras$n3
    }

    # if(is.unequal_clusters(paras$n2)) {
    #     paras$n3 <- length(unlist(paras$n2))
    #     paras_tx$n3 <- length(unlist(paras$n2))
    # }
    if(paras$partially_nested) {
        paras$sigma_cluster_intercept <- 0L
        paras$cor_cluster <- 0L
        paras$sigma_cluster_slope <- 0L
    }
    if(is.per_treatment(paras$dropout)) {
        paras_tx$dropout <- paras$dropout[[1]][[1]]
        paras$dropout <- paras$dropout[[1]][[2]]
    }

    list(control = paras,
         treatment = paras_tx)
}






#' Setup unbalanced cluster sizes
#'
#' Helps specifying unequal cluster sizes with \code{\link{study_parameters}}
#'
#' @param ... Any number of separate numeric arguments specifying
#' each cluster's size
#' @param func A function that generates cluster sizes, used instead of \code{...}. See \emph{Details}.
#' @param trunc Cutoff for values generated by \code{func}, \code{x < trunc} are replaced,
#' used to avoid negative or 0 values.
#' @param replace Indicates what value to replace cluster sizes less than \code{trunc} with.
#'
#' @details
#' If \code{func} is used together with a function that generates random draws, e.g.
#' \code{rnorm} or \code{rpois}, then cluster sizes  (and possible the number of clusters),
#' will be threated as a random variable. This is mostly intended for simulations. However, if a
#' random function is used with \code{get_power}, power will
#' vary each time the function is called, and you would have to average over repeated calls
#' to get the expected power. Functions that output decimal numbers will be rounded
#' to the closest integer.
#'
#' @return An object of type 'plcp_unequal_clusters'
#' @seealso \code{\link{per_treatment}}
#' @export
#'
#' @examples
#' library(dplyr)
#' n2 <- unequal_clusters(5, 10, 15, 40)
#' p <- study_parameters(n1 = 11,
#'                       n2 = n2,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       sigma_error = 1,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)
#'
#' # verify cluster sizes
#' d <- simulate_data(p)
#' d %>%
#'     filter(time == 0) %>%
#'     group_by(treatment, cluster) %>%
#'     summarise(n = n())
#'
#' # Poisson distributed cluster sizes
#' n2 <- unequal_clusters(func = rpois(n = 5, lambda = 5))
#' p <- study_parameters(n1 = 11,
#'                       n2 = n2,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       sigma_error = 1,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)
#'
#' # Use per_treatment() to specify per treatment ------------------------------
#' n2 <- per_treatment(unequal_clusters(2, 2, 2, 2, 3, 4, 5),
#'                      unequal_clusters(10, 15))
#' p <- study_parameters(n1 = 11,
#'                       n2 = n2,
#'                       n3 = 3,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)
#'
#' # verify cluster sizes
#' d <- simulate_data(p)
#' d %>%
#'     filter(time == 0) %>%
#'     group_by(treatment, cluster) %>%
#'     summarise(n = n())
unequal_clusters <- function(..., func = NULL, trunc = 1, replace = 1) {
    if(length(list(...)) > 0 & !is.null(func)) stop("Can't combine manual cluster sizes and 'func'.")
    if(is.null(func)) {
        x <- list(cluster_sizes = ...)
        tmp <- "manual"
    } else {
        x <- as.list(match.call())$func
        tmp <- x
    }

    out <- function() {
        if(is.call(x)) x <- eval(x)
        x <- unlist(x)
        attr(x, "func") <- tmp
        attr(x, "trunc") <- trunc
        attr(x, "replace") <- replace
        x
    }

    class(out) <- "plcp_unequal_clusters"
   # x <- list(unequal_clusters=x)
    #class(x) <- "unequal_clusters"
    list(out)
}
is.unequal_clusters <- function(x) {
    if(is.per_treatment(x)) {
        tx <- x[[1]]$treatment[[1]]
        cc <- x[[1]]$control[[1]]

        res <- any(c(class(tx), class(cc)) == "plcp_unequal_clusters")
    } else res <- class(x[[1]]) == "plcp_unequal_clusters"

    res
}

#' Setup parameters that differ per treatment group
#'
#' Helps specifying unequal cluster sizes with \code{\link{study_parameters}},
#' e.g. different number of clusters in the treatment and control arm, or
#' different dropout patterns.
#'
#' @param control Value used for control group
#' @param treatment Value used for treatment group
#'
#' @details The type of object passed to \code{control} and \code{treatment}
#' will depend on the parameters in \code{\link{study_parameters}} that should
#' have different values per treatment group.
#'
#' @return An object of class "plcp_per_treatment"
#' @seealso \code{\link{unequal_clusters}}, \code{\link{study_parameters}},
#' \code{\link{dropout_weibull}}
#' @export
#'
#' @examples
#' n2 <- per_treatment(control = 10,
#'                     treatment = 20)
#' p <- study_parameters(n1 = 11,
#'                       n2 = n2,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)

per_treatment  <- function(control, treatment) {
    dots <- list(treatment = treatment,
                 control = control)
    #res <- lapply(dots, list)

    class(dots) <- "plcp_per_treatment"
    x <- list(dots)
    #class(x) <- "plcp_per_treatment"

    x
}
is.per_treatment <- function(x) {
    if(!is.function(x)) {
        return(class(x[[1]]) == "plcp_per_treatment")
    } else return(FALSE)
}
as.plcp <- function(.p) {
    paras <- .p
    if(is.data.frame(paras)) {
        paras <- as.list(paras)
        paras <- do.call(study_parameters, paras)
        #class(paras) <- append(c("plcp"), class(paras))
    }
    paras
}

#' Update a \code{study_parameters}-object with new settings
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param ... Any number of named arguments that should be updated
#' @details Currently only the arguments used to construct the original object
#' can be updated.
#'
#' @examples
#' p <- study_parameters(n1 = 11,
#'                       n2 = 10,
#'                       n3 = 3,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       cohend = -0.8)
#'
#' p <- update(p, icc_slope = 0.1)
#' get_ICC_slope(p)
#'
#' \dontrun{
#' # Using a "new" argument does not work (yet)
#' update(p, sigma_cluster_slope = 2)
#' }
#' @method update plcp
#' @export
update.plcp <- function(object, ...) {
    paras <- object
    args <- attr(paras, "call")

    new_args <- list(...)
    new <- check_new_argument(args, new_args)
    if(length(new) > 0) stop(paste0("Updating new arguments is not yet implemented. '", new, "' was not used in original call."))
    for(i in seq_along(new_args)) {
        args[[names(new_args[i])]] <- new_args[[i]]
    }

    paras <- do.call(study_parameters, args)

    paras
}
#' @export
update.plcp_multi <- function(object, ...) {
    update.plcp(object, ...)
}

check_new_argument <- function(args, new) {
   x <- lapply(args[names(new)], is.null)
   x <- which(unlist(x))

   names(x)
}

get_n2 <- function(paras, n1 = 1) {
    UseMethod("get_n2")
}

get_n2.plcp <- function(paras) {
    tmp <- prepare_paras(paras)

    n2_cc <- unlist(tmp$control$n2)
    if(tmp$control$partially_nested) {
        if(length(n2_cc) == 1) {
            n2_cc<- tmp$control$n3 * n2_cc
            } else n2_cc <- sum(n2_cc)
    }
    n2_tx <- unlist(tmp$treatment$n2)

    list(treatment = n2_tx,
               control = n2_cc)
}

get_n2_ <- function(paras) {
    n2 <- unlist(paras$n2)
    if(length(n2) == 1) {
        n3 <- paras$n3
    } else {
        n3 <- length(n2)
    }

    n3
}


get_n3 <- function(paras, n1 = 1) {
    UseMethod("get_n3")
}

get_n3.plcp <- function(paras) {
    tmp <- prepare_paras(paras)

    n3_cc <- get_n3_(tmp$control)
    if(tmp$control$partially_nested) n3_cc <- 0L
    n3_tx <- get_n3_(tmp$treatment)

    data.frame(treatment = n3_tx,
               control = n3_cc,
               total = n3_tx + n3_cc)
}
get_n3_ <- function(paras) {
    n2 <- unlist(paras$n2)
    if(length(n2) == 1) {
        n3 <- paras$n3
    } else {
        n3 <- length(n2)
    }

    unlist(n3)
}

get_tot_n <- function(paras, n = 1) {
    UseMethod("get_tot_n")
}
get_tot_n.plcp <- function(paras, n = NULL) {
    tmp <- prepare_paras(paras)
    paras_cc <- tmp$control
    paras_tx <- tmp$treatment

    n_cc <- get_tot_n_(paras_cc)
    n_tx <- get_tot_n_(paras_tx)

    data.frame(treatment = n_tx,
               control = n_cc,
               total = n_tx + n_cc)
}
get_tot_n.plcp_multi <- function(paras, n = 1) {
    get_tot_n.plcp(as.plcp(paras[n, ]))
}
get_tot_n_ <- function(paras) {
    n2 <- unlist(paras$n2)
    if(is.unequal_clusters(n2)) n2 <- n2[[1]]()

    if(length(n2) == 1) {
        tot_n <- paras$n3 * n2
    } else {
        tot_n <- sum(n2)

    }
}
