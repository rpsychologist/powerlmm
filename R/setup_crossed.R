    
solve_subject_intercept <- function(sigma_subject_intercept, tmp) {
    if (is.null(sigma_subject_intercept)) {
        icc_cluster <- tmp$icc_pre_cluster
        icc_cluster[is.na(icc_cluster)] <- 0
        tmp$sigma_subject_intercept <-
            sqrt(
                tmp$sigma_error^2 /
                    (1 - (tmp$icc_pre_subject)) *
                    (tmp$icc_pre_subject - icc_cluster)
            )
    }
    tmp
}
solve_subject_slope <- function(sigma_subject_slope, tmp) {
    if (is.null(sigma_subject_slope)) {
        icc <- tmp$icc_slope
        icc[is.na(icc) | is.null(icc)] <- 0
        tmp$sigma_subject_slope <- sqrt(tmp$var_ratio *
            tmp$sigma_error^2 * (1 - icc))
    }
    tmp
}

#' @param crossed logical; indicates if design is three-level crossed. If true
#' also solves for sigma_cluster_intercept_crossed
solve_cluster_intercept <- function(
    sigma_subject_intercept,
    sigma_cluster_intercept,
    sigma_cluster_intercept_crossed = NULL,
    tmp, 
    crossed = FALSE) {
    if (is.null(sigma_cluster_intercept)) {
        # from icc_pre_subject
        if (is.null(sigma_subject_intercept)) {
            tmp$sigma_cluster_intercept <-
                sqrt(
                    tmp$sigma_error^2 /
                        (1 - (tmp$icc_pre_subject)) *
                        tmp$icc_pre_cluster
                )
            # from sigma_subject_intercept
        } else {
            v0 <- with(tmp, 
            (icc_pre_cluster * (sigma_error^2 + sigma_subject_intercept^2))
            / (1 - icc_pre_cluster))
            tmp$sigma_cluster_intercept <- sqrt(v0)
        }
    }
    if (crossed && is.null(sigma_cluster_intercept_crossed)) {
        icc_cluster_crossed <- tmp$icc_pre_cluster_crossed
        icc_cluster_crossed[is.na(icc_cluster_crossed) || is.null(icc_cluster_crossed)] <- 0    
        tot_cluster_var <- tmp$sigma_cluster_intercept^2
        tmp$sigma_cluster_intercept <-
            sqrt(tot_cluster_var * (1 - icc_cluster_crossed))
        tmp$sigma_cluster_intercept_crossed <-
            sqrt(tot_cluster_var * icc_cluster_crossed)
    }
    tmp
}

solve_cluster_slope <- function(
    icc_slope,
    sigma_subject_slope,
    sigma_cluster_slope,
    sigma_cluster_slope_crossed = NULL,
    tmp,
    crossed = FALSE) {
    if (is.null(sigma_cluster_slope)) {
        # check if NA
        if (is.null(icc_slope) || all(is.na(icc_slope))) {
            tmp$sigma_cluster_slope <- NA
        } else {
            # solve from icc_slope and var_ratio
            if (is.null(sigma_subject_slope)) {
                v1 <- with(tmp, var_ratio * sigma_error^2 * icc_slope)
                tmp$sigma_cluster_slope <- sqrt(v1)
                # solve from subject_slope and icc_slope
            } else {
                x <- with(tmp, sigma_subject_slope^2 / (1 - icc_slope))
                v1 <- x - tmp$sigma_subject_slope^2
                v1 <- sqrt(v1)
                tmp$sigma_cluster_slope <- v1
            }
        }
    }
    if (crossed && is.null(sigma_cluster_slope_crossed)) {
        icc_cluster_crossed <- tmp$icc_slope_crossed
        icc_cluster_crossed[is.na(icc_cluster_crossed) || is.null(icc_cluster_crossed)] <- 0
        tot_cluster_var <- tmp$sigma_cluster_slope^2
        tmp$sigma_cluster_slope <-
            sqrt(tot_cluster_var * (1 - icc_cluster_crossed))
        tmp$sigma_cluster_slope_crossed <-
            sqrt(tot_cluster_var * tmp$icc_slope_crossed)
    }
    tmp
}


#' Setup crossed study parameters
#'
#' Setup the parameters for calculating power for longitudinal multilevel studies
#' comparing two groups. Ordinary two-level models (subjects with repeated measures),
#' and longitudinal three-level models with clustering due to therapists, schools, provider etc,
#' are supported. Random slopes at the subject level and cluster level are
#' possible. Cluster sizes can be unbalanced, and vary by treatment.
#' Partially nested designs are supported. Missing data can also be accounted
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
#' @param icc_slope_crossed Proportion of \code{icc_slope} slope variance
#' at the cluster level casued by the crossed cluster:treatment interaction.
#' @param var_ratio Ratio of the random
#' slope variance to the within-subject variance.
#' @param icc_pre_subject Proportion of baseline
#' variance at the subject level. N.B. the variance at the subject-level also
#' included the cluster-level variance. If there's no random slopes, this would
#' be the subject-level ICC, i.e. correlation between time points.
#' @param icc_pre_cluster Proportion of baseline
#' variance at the cluster level.
#' @param icc_pre_cluster_crossed Proportion of \code{icc_pre_cluster} variance related to the
#' cluster:treatment interaction.
#' @param cor_subject Correlation between the subject-level random intercept
#'  and slopes.
#' @param cor_cluster A 3 x 3 matrix of the cluster-level variance components.
#' @param cor_within Correlation of the level 1 residual. Currently ignored in
#' the analytical power calculations.
#' @param effect_size The treatment effect. Either a \code{numeric} indicating the mean
#' difference (unstandardized) between the treatments at posttest, or a standardized effect
#' using the \code{\link{cohend}} helper function.
#' @param cohend \emph{Deprecated}; now act as a shortcut to \code{\link{cohend}} helper function.
#' Equivalent to using \code{effect_size = cohend(cohend, standardizer = "pretest_SD", treatment = "control")}
#' @param partially_nested \code{logical}; indicates if there's clustering in both
#' arms or only in the treatment arm.
#' @param dropout Dropout process, see \code{\link{dropout_weibull}} or
#' \code{\link{dropout_manual}}. Assumed to be 0 if \code{NULL}.
#' @param deterministic_dropout \code{logical}; if \code{FALSE} the input to
#' \code{dropout} will be treated as random and dropout will be sampled
#' from a multinomial distribution. N.B.: the random dropout will be
#' sampled independently in both treatment arms.
#' @return A \code{list} or \code{data.frame} of parameters values, either of
#' class \code{plcp} or \code{plcp_multi} if multiple parameters are compared.
#'
#' @inherit study_parameters details
#'
#' @details
#' \bold{Crossed Factor}
#' The interpretation of the crossed variance components depends on the coding of the
#' \code{treatment} variable.
#' 
#' @seealso \code{\link{cohend}}, \code{\link{get_power}}, \code{\link{simulate.plcp}}
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
#'                       effect_size = cohend(-0.8))
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
#'                       effect_size = cohend(-0.8))
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
#'                       effect_size = cohend(-0.8))
#'
#' get_power(p)
#'
#' ## Two-level model with subject-level random slope
#' p <- study_parameters(n1 = 11,
#'                       n2 = 40,
#'                       icc_pre_subject = 0.5,
#'                       var_ratio = 0.03,
#'                       effect_size = cohend(-0.8))
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
#'                       effect_size = cohend(c(-0.5, -0.8))
#'                       )
#'
#' get_power(p)
#' @export
study_parameters.plcp_design_crossed <- function(design,
                                     n1,
                                     n2,
                                     n3 = 1,
                                     T_end = NULL,
                                     fixed_intercept = 0L,
                                     fixed_tx = 0L,
                                     fixed_slope = 0L,
                                     sigma_subject_intercept = NULL,
                                     sigma_subject_slope = NULL,
                                     sigma_cluster_intercept = NULL,
                                     sigma_cluster_slope = NULL,
                                     sigma_cluster_intercept_crossed = NULL,
                                     sigma_cluster_slope_crossed = NULL,
                                     sigma_error = 10,
                                     cor_subject = 0L,
                                     cor_cluster_intercept_slope = 0L,
                                     cor_cluster_intercept_intercept_tx = 0L,
                                     cor_cluster_intercept_slope_tx = 0L,
                                     cor_cluster_slope_intercept_tx = 0L,
                                     cor_cluster_slope_slope_tx = 0L,
                                     cor_cluster_intercept_tx_slope_tx = 0L,
                                     cor_within = 0L,
                                     var_ratio = NULL,
                                     icc_slope = NULL,
                                     icc_slope_crossed = NULL,
                                     icc_pre_subject = NULL,
                                     icc_pre_cluster = NULL,
                                     icc_pre_cluster_crossed = NULL,
                                     effect_size = 0L,
                                     cohend = NULL,
                                     dropout = 0L,
                                     deterministic_dropout = TRUE,
                                     ...) {
    print("crossed")
    # deprecated Cohen's d
    if(!is.null(cohend)) {
        effect_size <- cohend(cohend, standardizer = "pretest_SD", treatment = "control")
    }
    if(is.per_treatment(n3)) stop("the number of clusters can't differ per preatment in a crossed design.", call. = FALSE)
    args <- list(
        n1 = n1,
        n2 = n2,
        n3 = n3,
        T_end = T_end,
        fixed_intercept = fixed_intercept,
        fixed_tx = fixed_tx,
        fixed_slope = fixed_slope,
        sigma_subject_intercept = sigma_subject_intercept,
        sigma_subject_slope = sigma_subject_slope,
        sigma_cluster_intercept = sigma_cluster_intercept,
        sigma_cluster_slope = sigma_cluster_slope,
        sigma_cluster_intercept_crossed = sigma_cluster_intercept_crossed,
        sigma_cluster_slope_crossed = sigma_cluster_slope_crossed,
        sigma_error = sigma_error,
        cor_subject = cor_subject,
        cor_cluster_intercept_slope = cor_cluster_intercept_slope,
        cor_cluster_intercept_intercept_tx = cor_cluster_intercept_intercept_tx,
        cor_cluster_intercept_slope_tx = cor_cluster_intercept_slope_tx,
        cor_cluster_slope_intercept_tx = cor_cluster_slope_intercept_tx,
        cor_cluster_slope_slope_tx = cor_cluster_slope_slope_tx,
        cor_cluster_intercept_tx_slope_tx = cor_cluster_intercept_tx_slope_tx,
        cor_within = cor_within,
        var_ratio = var_ratio,
        icc_slope = icc_slope,
        icc_slope_crossed = icc_slope_crossed,
        icc_pre_subject = icc_pre_subject,
        icc_pre_cluster = icc_pre_cluster,
        icc_pre_cluster_crossed = icc_pre_cluster_crossed,
        effect_size = effect_size,
        dropout = dropout,
        deterministic_dropout = deterministic_dropout,
        family = design$family
    )
    save_call <- args
    save_call$design <- design
    tmp_args <- args[!vapply(args, is.null, logical(1))]
    ## Default NA
    if(is.null(icc_pre_subject) & is.null(sigma_subject_intercept)) {
        tmp_args$icc_pre_subject <- NA
    }
    if (is.null(var_ratio) & (is.null(sigma_subject_slope) || is.na(sigma_subject_slope)) & is.null(sigma_cluster_slope)) {
        tmp_args$var_ratio <- NA
    }
    if(is.null(icc_pre_cluster) & is.null(sigma_cluster_intercept)) {
        tmp_args$icc_pre_cluster <- NA
    }
    if(is.null(icc_pre_cluster_crossed) & is.null(sigma_cluster_intercept_crossed)) {
        tmp_args$icc_pre_cluster_crossed <- NA
    }
    if(is.null(icc_slope) & is.null(sigma_cluster_slope)) {
        tmp_args$icc_slope <- NA
    }
    if(is.null(icc_slope_crossed) & is.null(sigma_cluster_slope_crossed)) {
        tmp_args$icc_slope_crossed <- NA
    }
     tmp <- expand.grid(tmp_args)
    # ## --
    #
    ## icc_cluster > icc_subject
    if(!is.null(icc_pre_cluster) & !is.null(icc_pre_subject)) {
        if(any((tmp$icc_pre_cluster > tmp$icc_pre_subject), na.rm=TRUE)) {
            stop("'icc_pre_cluster' can't be larger than 'icc_pre_subject'", call. = FALSE)
        }
    }
    #
    # # check cluster slope variance exists when var ratio is NA or zero.
    # if(!is.null(tmp$var_ratio)) {
    #     if(is.na(tmp$var_ratio) && any(tmp$icc_slope >= 0, na.rm=TRUE)) {
    #         stop("'icc_slope' can't be >= 0 when 'var_ratio' or 'sigma_subject_slope' is NA", call. = FALSE)
    #     }
    #     if((any(tmp$var_ratio == 0, na.rm=TRUE) | any(sigma_subject_slope == 0)) &&
    #        any(tmp$icc_slope > 0, na.rm=TRUE)) {
    #         stop("'icc_slope' can't be > 0 when 'var_ratio' or 'sigma_subject_slope' is 0", call. = FALSE)
    #     }
    # }

    ## Default T_end
    if(is.null(args$T_end)) tmp$T_end <- tmp$n1 - 1


    ## Solve raw values
    # Solve subject_intercept
    tmp <- solve_subject_intercept(sigma_subject_intercept, tmp)
    # solve subject_slope
    tmp <- solve_subject_slope(sigma_subject_slope, tmp)
    # Solve cluster_intercept
    tmp <- solve_cluster_intercept(
        sigma_subject_intercept = sigma_subject_intercept,
        sigma_cluster_intercept = sigma_cluster_intercept,
        sigma_cluster_intercept_crossed = sigma_cluster_intercept_crossed,
        tmp = tmp,
        crossed = TRUE
    )
    # Solve cluster_slope
    tmp <- solve_cluster_slope(
        icc_slope = icc_slope,
        sigma_subject_slope = sigma_subject_slope,
        sigma_cluster_slope = sigma_cluster_slope,
        sigma_cluster_slope_crossed = sigma_cluster_slope_crossed,
        tmp = tmp,
        crossed = TRUE
    )
    # keep cols
    cols <- which(colnames(tmp) %in% c(
        "icc_slope",
        "icc_slope_crossed",
        "var_ratio",
        "icc_pre_cluster",
        "icc_pre_cluster_crossed",
        "icc_pre_subject")
    )
    cols <- colnames(tmp)[cols]
    paras <- tmp[, !(names(tmp) %in% cols)]


    # Single or multi?
    paras <- .make_single_or_multi(paras, model = "crossed")

    paras$design <- "plcp_design_crossed"
    # avoid problems with prepare_paras()
    paras$partially_nested <- FALSE
    paras$custom_model <- FALSE
    # if(is.data.frame(paras)) {
    #     class(paras) <- append(c("plcp_multi"), class(paras))
    # } else class(paras) <- append(c("plcp_crossed", "plcp"), class(paras))
    # #
    # # Default cor_*
    # if(is.null(paras$cor_cluster)) paras$cor_cluster <- cor_cluster
    # if(is.null(paras$cor_subject)) paras$cor_subject <- cor_subject

    # # Classes
     class(paras) <- append(class(paras), c("plcp_crossed"))
    # if(all(is.na(paras$sigma_cluster_slope)) &
    #    all(is.na(paras$sigma_cluster_intercept))) {
    #     class(paras) <- append(class(paras), c("plcp_2lvl"))
    # } else if (all(!is.na(paras$sigma_cluster_slope)) |
    #            all(!is.na(paras$sigma_cluster_intercept))) {
    #     class(paras) <- append(class(paras), c("plcp_3lvl"))
    # } else {
    #     class(paras) <- append(class(paras), c("plcp_mixed"))
    # }

    attr(paras, "call") <- save_call
    paras

}

