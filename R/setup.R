#' Setup study parameters
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
#' @param var_ratio Ratio of the random
#' slope variance to the within-subject variance.
#' @param icc_pre_subject Amount of baseline
#' variance at the subject level. N.B. the variance at the subject-level also
#' included the cluster-level variance. If there's no random slopes, this would
#' be the subject-level ICC, i.e. correlation between time points.
#' @param icc_pre_cluster Amount of baseline
#' variance at the cluster level.
#' @param cor_subject Correlation between the subject-level random intercept
#'  and slopes.
#' @param cor_cluster Correlation between the cluster-level random intercept
#' and slopes.
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
#' solved for. Only the cluster-level parameters can be solved when standardized and
#' raw values are mixed. \code{sigma_error} is 10 by default. More information regarding
#' the standardized inputs are available in the two-level and three-level vignettes.
#'
#' \bold{Difference between 0 and NA}
#'
#' For the variance components \code{0} and \code{NA/NULL} have different meanings.
#' A parameter that is 0 is still kept in the model, e.g. if \code{icc_pre_cluster = 0}
#' a random intercept is estimated at the cluster level, but the true value is 0.
#' If the argument is either \code{NULL} or \code{NA} it is excluded from the model.
#' This choice will matter when running simulations, or if Satterthwaite \emph{dfs} are used.
#'
#' The default behavior if a parameters is not specified is that \code{cor_subject} and
#' \code{cor_cluster} is 0, and the other variance components are \code{NULL}.
#'
#' \bold{Effect size and Cohen's d}
#'
#' The argument \code{effect_size} let's you specify the average difference in change
#' between the treatment groups. You can either pass a \code{numeric} value to define
#' the raw difference in means at posttest, or use a standardized effect size, see
#' \code{\link{cohend}} for more details on the standardized effects.
#'
#' The argument \code{cohend} is kept for legacy reasons, and is equivalent to using
#' \code{effect_size = cohend(cohend, standardizer = "pretest_SD", treatment = "control")}.
#'
#' \bold{Two- or three-level models}
#'
#' If either \code{sigma_cluster_slope} or \code{icc_slope} and
#'  \code{sigma_cluster_intercept} or \code{icc_pre_cluster} is
#' \code{NULL} it will be assumed a two-level design is wanted.
#'
#' \bold{Unequal cluster sizes and unbalanced allocation}
#'
#' It is possible to specify different cluster sizes using
#' \code{\link{unequal_clusters}}. Cluster sizes can vary between treatment arms
#' by also using \code{\link{per_treatment}}. The number of clusters per treatment can
#' also be set by using \code{\link{per_treatment}}. Moreover, cluster
#' sizes can be sampled from a distribution, and treated as a random variable.
#' See \code{\link{per_treatment}} and \code{\link{unequal_clusters}} for examples of their use.
#'
#' \bold{Missing data and dropout}
#'
#' Accounting for missing data in the power calculations is possible. Currently,
#' \code{dropout} can be specified using either \code{\link{dropout_weibull}} or
#' \code{\link{dropout_manual}}. It is possible to have different dropout
#' patterns per treatment group using \code{\link{per_treatment}}. See their
#' respective help pages for examples of their use.
#'
#' If \code{deterministic_dropout = TRUE} then the proportion of dropout is treated is fixed.
#' However, exactly which subjects dropout is randomly sampled within treatments. Thus,
#' clusters can become slightly unbalanced, but generally power varies little over realizations.
#'
#' For \emph{random dropout}, \code{deterministic_dropout = FALSE}, the proportion
#' of dropout is converted to the probability of having exactly \emph{i} measurements,
#' and the actual dropout is sampled from a multinomial distribution. In this case, the proportion of
#' dropout varies over the realizations from the multinomial distribution, but will
#' match the dropout proportions in expectation. The random dropout in
#' each treatment group is sampled from independent multinomial distributions.
#'
#' Generally, power based on fixed dropout is a good approximation of random dropout.
#'
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
study_parameters <- function(n1,
                             n2,
                             n3 = 1,
                             T_end = NULL,
                             fixed_intercept = 0L,
                             fixed_slope = 0L,
                             sigma_subject_intercept = NULL,
                             sigma_subject_slope = NULL,
                             sigma_cluster_intercept = NULL,
                             sigma_cluster_slope = NULL,
                             sigma_error = 10,
                             cor_subject = 0L,
                             cor_cluster = 0L,
                             cor_within = 0L,
                             var_ratio = NULL,
                             icc_slope = NULL,
                             icc_pre_subject = NULL,
                             icc_pre_cluster = NULL,
                             effect_size = 0L,
                             cohend = NULL,
                             partially_nested = FALSE,
                             dropout = 0L,
                             deterministic_dropout = TRUE) {

    #if(!is.per_treatment(n2) & length(n2) == 1) n2 <- list(n2)

    # deprecated Cohen's d
    if(!is.null(cohend)) {
        effect_size <- cohend(cohend, standardizer = "pretest_SD", treatment = "control")
    }

    # drop out checks
    if(is.numeric(dropout) && any(dropout != 0)) stop("'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'", call. = FALSE)
    if(is.per_treatment(dropout)) {
        tx <- dropout[[1]]$treatment
        cc <- dropout[[1]]$control
        if(is.numeric(cc) && any(cc != 0)) stop("Control group's 'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'", call. = FALSE)
        if(is.numeric(tx) && any(tx != 0)) stop("Treatment group's 'dropout' should be 0 or created by 'dropout_manual' or 'dropout_weibull'", call. = FALSE)
    }

    # warn n3 is ignored
    if(is.unequal_clusters(n2) & !is.per_treatment(n2) & is.per_treatment(n3)) {
        message("'n3' per_treatment argument is ignored. 'n3' is automatically based on length of unequal_clusters")
    }

    # icc_slope + sigma_*_slope
    if(!is.null(icc_slope) & (!is.null(sigma_subject_slope) &
                              !is.null(sigma_cluster_slope))) {
        stop("Can't use 'icc_slope' with both 'sigma_subject_slope' and 'sigma_cluster_slope'", call. = FALSE)

    }
    if(!is.null(sigma_subject_slope) &
       !is.null(var_ratio)) {
        stop("Can't use both 'sigma_subject_slope' and 'var_ratio'", call. = FALSE)

    }
    if(!is.null(sigma_cluster_slope) &
       !is.null(var_ratio)) {
        stop("Can't use both 'sigma_cluster_slope' and 'var_ratio'", call. = FALSE)

    }
    if(!is.null(sigma_subject_intercept) & !is.null(icc_pre_subject)) {
        stop("Can't use both 'icc_pre_subject' and 'sigma_subject_intercept'", call. = FALSE)
    }

    if(!is.null(sigma_cluster_intercept) & !is.null(icc_pre_cluster)) {
        stop("Can't use both 'icc_pre_cluster' and 'sigma_cluster_intercept'", call. = FALSE)
    }

    if(is.null(sigma_subject_intercept) & is.null(icc_pre_subject)) {
           stop("Both 'sigma_subject_intercept' and 'icc_pre_subject' can't be NULL", call. = FALSE)
    }
    if(!is.null(sigma_cluster_slope) & !is.null(icc_slope)) {
        stop("Can't use both 'icc_slope' and 'sigma_cluster_slope'", call. = FALSE)
    }
    if(!is.null(sigma_cluster_slope) & is.null(sigma_subject_slope)) {
        stop("'sigma_subject_slope' is missing", call. = FALSE)
    }
    if(is.null(icc_pre_cluster) &
       !is.null(sigma_cluster_intercept) &
       !is.null(icc_pre_subject)) {
        stop("'icc_pre_subject' and 'sigma_cluster_intercept' can't be combined, use 'icc_pre_cluster'", call. = FALSE)
    }


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
        effect_size = effect_size,
        partially_nested = partially_nested,
        dropout = dropout,
        deterministic_dropout = deterministic_dropout
    )
    save_call <- args

    tmp_args <- args[!vapply(args, is.null, logical(1))]


    ## Default NA
    if(is.null(icc_pre_subject) & is.null(sigma_subject_intercept)) {
        tmp_args$icc_pre_subject <- NA
    }
    if(is.null(var_ratio) & (is.null(sigma_subject_slope) || is.na(sigma_subject_slope)) & is.null(sigma_cluster_slope)) {
        tmp_args$var_ratio <- NA
    }
    if(is.null(icc_pre_cluster) & is.null(sigma_cluster_intercept)) {
        tmp_args$icc_pre_cluster <- NA
    }
    if(is.null(icc_slope) & is.null(sigma_cluster_slope)) {
        tmp_args$icc_slope <- NA
    }
    tmp <- expand.grid(tmp_args)
    ## --

    ## icc_cluster > icc_subject
    if(!is.null(icc_pre_cluster) & !is.null(icc_pre_subject)) {
        if(any((tmp$icc_pre_cluster > tmp$icc_pre_subject), na.rm=TRUE)) {
            stop("'icc_pre_cluster' can't be larger than 'icc_pre_subject'", call. = FALSE)
        }
    }

    # check cluster slope variance exists when var ratio is NA or zero.
    if(!is.null(tmp$var_ratio)) {
        if(is.na(tmp$var_ratio) && any(tmp$icc_slope >= 0, na.rm=TRUE)) {
            stop("'icc_slope' can't be >= 0 when 'var_ratio' or 'sigma_subject_slope' is NA", call. = FALSE)
        }
        if((any(tmp$var_ratio == 0, na.rm=TRUE) | any(sigma_subject_slope == 0)) &&
           any(tmp$icc_slope > 0, na.rm=TRUE)) {
            stop("'icc_slope' can't be > 0 when 'var_ratio' or 'sigma_subject_slope' is 0", call. = FALSE)
        }
    }

    ## Default T_end
    if(is.null(args$T_end)) tmp$T_end <- tmp$n1 - 1


    ## Solve raw values

    # Solve subject_intercept
    if(is.null(sigma_subject_intercept)) {

        icc_cluster <- tmp$icc_pre_cluster
        icc_cluster[is.na(icc_cluster)] <- 0

        tmp$sigma_subject_intercept <-
            sqrt(
                tmp$sigma_error^2 /
                    (1-(tmp$icc_pre_subject)) *
                    (tmp$icc_pre_subject - icc_cluster)
            )
    }

    # solve subject_slope
    if(is.null(sigma_subject_slope)) {

        icc <- tmp$icc_slope
        icc[is.na(icc) | is.null(icc)] <- 0
        tmp$sigma_subject_slope <- sqrt(tmp$var_ratio *
                                                 tmp$sigma_error^2 * (1-icc))

    }
    # Solve cluster_intercept
    if(is.null(sigma_cluster_intercept)) {
        # from icc_pre_subject
        if(is.null(sigma_subject_intercept)) {
            tmp$sigma_cluster_intercept <-
                sqrt(
                    tmp$sigma_error^2 /
                        (1-(tmp$icc_pre_subject)) *
                        tmp$icc_pre_cluster
                )
        # from sigma_subject_intercept
        } else {
            v0 <- with(tmp, (icc_pre_cluster*(sigma_error^2 + sigma_subject_intercept^2))
                       /(1-icc_pre_cluster))
            tmp$sigma_cluster_intercept <- sqrt(v0)
        }

    }

    # Solve cluster_slope
    if(is.null(sigma_cluster_slope)) {

        # check if NA
        if(is.null(icc_slope) || all(is.na(icc_slope)) ) {

            tmp$sigma_cluster_slope <- NA
        } else {
            # solve from icc_slope and var_ratio
            if(is.null(sigma_subject_slope)) {
                v1 <- with(tmp, var_ratio * sigma_error^2 * icc_slope)
                tmp$sigma_cluster_slope <- sqrt(v1)
                # solve from subject_slope and icc_slope
            } else {
                x <- with(tmp, sigma_subject_slope^2/(1-icc_slope))
                v1 <- x - tmp$sigma_subject_slope^2
                v1 <- sqrt(v1)
                tmp$sigma_cluster_slope <- v1

            }
        }
    }


    # keep cols
    cols <- which(colnames(tmp) %in% c("icc_slope",
                                       "var_ratio",
                                       "icc_pre_cluster",
                                       "icc_pre_subject"))
    cols <- colnames(tmp)[cols]
    paras <- tmp[, !(names(tmp) %in% cols)]


    # Single or multi?
    if((is.data.frame(paras) & nrow(paras) == 1)) {
        paras <- as.list(paras)
    }
    if(is.data.frame(paras)) {
        class(paras) <- append(c("plcp_multi"), class(paras))
    } else class(paras) <- append(c("plcp"), class(paras))

    # Default cor_*
    if(is.null(paras$cor_cluster)) paras$cor_cluster <- cor_cluster
    if(is.null(paras$cor_subject)) paras$cor_subject <- cor_subject

    # Classes
    if(all(is.na(paras$sigma_cluster_slope)) &
       all(is.na(paras$sigma_cluster_intercept))) {
        class(paras) <- append(class(paras), c("plcp_2lvl"))
    } else if (all(!is.na(paras$sigma_cluster_slope)) |
               all(!is.na(paras$sigma_cluster_intercept))) {
        class(paras) <- append(class(paras), c("plcp_3lvl"))
    } else {
        class(paras) <- append(class(paras), c("plcp_mixed"))
    }

    attr(paras, "call") <- save_call
    paras

}

sim_parameters <- function(...) {
    dots <- list(...)
    warning("sim_parameters is deprecated")
    do.call(study_parameters, dots)
}


print_per_treatment <- function(n, width = 0, n2 = FALSE, hanging = 19) {

    x <- lapply(seq_along(n), function(i) {
        print_per_treatment_(i, x = n, n2 = n2)$lab
    })
    x <- format(x, width = width)
    x <- paste(x, " (", names(n), ")", sep ="")
    collapse <- paste0("\n", paste(rep(" ", hanging), collapse = ""))
    x <- paste(unlist(x), collapse = collapse)
    x
}
print_per_treatment_ <- function(i, x, n2 = FALSE) {
    name <- names(x)[i]
    x <- x[[i]]
    if(n2 & length(unique(x)) == 1) {
        if(attr(x, "func")) {
            x_num <- NA
            x_lab <- paste(unique(x))
        } else {
            x_num <- unique(x)
            if(length(x) > 1) {
                x_lab <- paste(unique(x),"x", length(x))
            } else {
                x_lab <- paste(unique(x))
            }

        }
    } else if(length(unique(x)) == 1) {
        x_lab <- x
        x_num <- x
    } else {
        x_lab <- x
        x_num <- NA
    }
    list("lab" = paste(paste(unlist(x_lab), collapse = ", "), sep =""),
         "num" = x_num)
}

deparse_n2 <- function(n2) {
    n2_attr <- attr(n2, "func")
    if(!is.null(n2_attr) && (n2_attr != "manual")) {
        n2 <- deparse(n2_attr)
        attr(n2, "func") <- TRUE
    } else attr(n2, "func") <- FALSE

    n2

}
truncate_n2 <- function(n2) {
    n2 <- gsub("^.*::", "", n2)
    n2 <- strtrim(n2, 20)

}
prepare_print_n2 <- function(x) {
    n2 <- get_n2(x)
    n2$treatment <- deparse_n2(n2$treatment)
    n2_attr <- attr(n2$control, "func")
    if(attr(n2$control, "per_treatment")) {
        n2$control <- deparse_n2(n2$control)

    } else if(!is.null(n2_attr) && (n2_attr != "manual")) {
        n2$control <- "-"
        attr(n2$control, "func") <- FALSE
    } else {
        n2$control <- deparse_n2(n2$control)
    }


    n2
}
prepare_print_plcp <- function(x, two_level = FALSE) {
    n1 <- x$n1
    n2 <- prepare_print_n2(x)
    n3 <- get_n3(x)
    tot_n <- get_tot_n(x)
    width <- max(nchar(print_per_treatment_(1, n2, n2 = TRUE)$lab),
                 nchar(print_per_treatment_(2, n2, n2 = TRUE)$lab),
                 nchar(print_per_treatment_(3, tot_n)$lab))
    if(two_level) width <- max(vapply(tot_n, nchar, numeric(1)))
    n2 <- print_per_treatment(n2, width = width, n2 = TRUE)
    n3 <- print_per_treatment(n3, width = width)

    tot_n <- print_per_treatment(tot_n, width = width)

    icc_slope <- round(get_ICC_slope(x), 2)
    var_ratio <- round(get_var_ratio(x), 2)
    icc_pre_clusters <- round(get_ICC_pre_clusters(x), 2)
    icc_pre_subjects <- round(get_ICC_pre_subjects(x), 2)

    effect <- get_effect_size(x)
    effect_label <- ifelse(effect$standardizer == "raw", "raw", "Cohen's d")
    if(x$partially_nested) {
        effect_label <- ifelse(effect_label == "Cohen's d",
                               paste(effect_label, " [SD: ", effect$standardizer,
                                     ", ", effect$treatment, "]", sep = ""),
                               effect_label)
    } else {
        effect_label <- ifelse(effect_label == "Cohen's d",
                               paste(effect_label, " [SD: ", effect$standardizer, "]", sep = ""),
                               effect_label)
    }

    effect <- paste(effect$ES, " (", effect_label,")", sep = "")

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
                          `effect_size` = effect,
                          method = "Study setup (three-level)"),
                     class = "power.htest")
    attr(res, "width") <- width
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

        if(!is.null(n2_attr) &
           n2_attr != "manual") res$note <- "n2 is randomly sampled"
        if(!is.null(n2_attr) &
           n2_attr != "manual" &
           !is.per_treatment(x$n2)) attr(res, "same_dist") <- TRUE
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
    if(!is.null(attr(res, "same_dist"))) message("The same random 'n2' sample is used for both treatment groups,\n'per_treatment' is needed to sample both groups independently.")

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


#' Return the raw difference between the groups at posttest
#'
#' Used internally to calculate the difference in change over time
#' between the two treatment groups.
#'
#' @param object A \code{\link{study_parameters}}-object.
#'
#' @return A \code{numeric} indicating the mean difference between the treatment and
#' control group at posttest.
#' @export
get_slope_diff <- function(object) {
    UseMethod("get_slope_diff")
}

#' @rdname get_slope_diff
#' @export
get_slope_diff.plcp <- function(object) {
    object$sigma_subject_intercept[is.na(object$sigma_subject_intercept)] <- 0
    object$sigma_cluster_intercept[is.na(object$sigma_cluster_intercept)] <- 0

    if(inherits(object$effect_size[[1]], "plcp_cohend")) {
        slope <- object$effect_size[[1]]$set(object)
    } else if(is.numeric(unlist(object$effect_size))) {
        slope <- unlist(object$effect_size)
    }

    slope
}
#' @rdname get_slope_diff
#' @export
get_slope_diff.plcp_multi <- function(object) {
    vapply(1:nrow(object), function(i) {
        p <- as.plcp(object[i, ])
        get_slope_diff.plcp(p)
    },
    numeric(1))

}



#' Use Cohen's d as the effect size in \code{study_parameters}
#'
#' This function is used as input to the \code{effect_size} argument in \code{study_parameters},
#' if standardized effect sizes should be used. The choice of the denominator differs between fields,
#' and this function supports the common ones: pre- or posttest SD, or the random slope SD.
#'
#' @param ES \code{numeric}; value of the standardized effect size. Can be a vector.
#' @param standardizer \code{character}; the standardizer (denominator) used to calculate
#' Cohen's d. Allows options are: "pretest_SD", "posttest_SD", or "slope_SD".
#' See Details from more information.
#' @param treatment \code{character}; indicates if the \code{standardizer} should
#' be based on the "treatment" or "control" group---this only matters for 3-level partially
#' nested designs.
#'
#' @details
#'
#' \strong{Standardizing using the \code{pretest_SD} or \code{posttest_SD}}
#'
#' For these effect sizes, ES indicates the standardized difference between
#' the treatment groups at posttest (\code{T_end}), standardized by using
#' either the implied standard deviation at pretest or posttest. Thus, the actual
#' raw differences in average slopes between the treatments are,
#'
#' \code{slope_diff = (ES * SD)/T_end}.
#'
#' \strong{\code{slope_SD}: standardizing using the random slopes}
#'
#' This standardization is quite different from using the pretest or posttest SD.
#' Here the average slope difference is standardized using the total SD of the random slopes.
#' This is done by e.g. Raudenbush and Liu (2001). \strong{NB}, for this effect size
#' \code{ES} indicates the difference in change per unit time, and not at posttest. Thus, the raw
#' difference in average slopes is,
#'
#' \code{slope_diff = ES * slope_SD}.
#'
#' For a 3-level model, \code{slope_SD = sqrt(sigma_subject_slope^2 + sigma_cluster_slope^2)}.
#'
#' @seealso \code{\link{study_parameters}}
#'
#' @references Raudenbush, S. W., & Liu, X. F. (2001). Effects of study duration,
#' frequency of observation, and sample size on power in studies of group differences
#' in polynomial change. \emph{Psychological methods}, 6(4), 387.
#'
#' @return A \code{list} of the same length as \code{ES}. Each element is a named list
#' of class \code{plcp_cohend}, with the elements:
#' \itemize{
#' \item \code{set}: A helper \code{function} that converts the standardized ES to raw values.
#' Accepts a \code{study_parameters} objects, and returns a \code{numeric} indicating the
#' raw difference between the treatment at posttest.
#' \item \code{get}: contains a list with the original call: "ES", "standardizer", and "treatment".
#' }
#'
#'
#' @export
#'
#' @examples
#'
#' # Pretest SD
#' p <- study_parameters(n1 = 11,
#'                       n2 = 20,
#'                       icc_pre_subject = 0.5,
#'                       cor_subject = -0.4,
#'                       var_ratio = 0.03,
#'                       effect_size = cohend(0.4, standardizer = "pretest_SD"))
#'
#' get_slope_diff(p)
#'
#' # using posttest SD,
#' # due to random slope SD will be larger at posttest
#' # thus ES = 0.4 indicate larger raw slope difference
#' # using posttest SD
#' p <- update(p, effect_size = cohend(0.4,
#'                                     standardizer = "posttest_SD"))
#' get_slope_diff(p)
#'
#'
#' # Random slope SD
#' p <- study_parameters(n1 = 11,
#'                       n2 = 20,
#'                       icc_pre_subject = 0.5,
#'                       cor_subject = -0.4,
#'                       var_ratio = 0.03,
#'                       effect_size = cohend(0.4, standardizer = "slope_SD"))
#'
#' # Partially nested ----------------------------------------------------------
#' p <- study_parameters(n1 = 11,
#'                       n2 = 20,
#'                       n3 = 4,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0.25,
#'                       cor_subject = -0.4,
#'                       var_ratio = 0.03,
#'                       partially_nested = TRUE,
#'                       effect_size = cohend(0.4, standardizer = "pretest_SD")
#'                       )
#' # Default is to use control groups SD
#' get_slope_diff(p)
#'
#' # Treatment group's SD also include cluster-level intercept variance.
#' # Thus, ES of 0.4 will indicate a larger raw difference
#' # using the treatment group's SD
#' p <- update(p, effect_size = cohend(0.4,
#'                                     standardizer = "pretest_SD",
#'                                     treatment = "treatment"))
#' get_slope_diff(p)
#'
#' ## Combine multiple values, and raw and standardized effects ----------------
#' p <- study_parameters(n1 = 11,
#'                       n2 = 20,
#'                       icc_pre_subject = 0.5,
#'                       cor_subject = -0.4,
#'                       var_ratio = 0.03,
#'                       effect_size = c(-5, 9,
#'                                       cohend(c(0.5, 0.8), standardizer = "pretest_SD"),
#'                                       cohend(c(0.5, 0.8), standardizer = "posttest_SD")))
#'
#'
#' ## Recreate results in Raudenbush & Liu 2001 --------------------------------
#' rauden_liu <- function(D, f, n = 238) {
#'     n1 <- f * D + 1
#'     p <- study_parameters(n1 = n1,
#'                           n2 = n/2,
#'                           T_end = D,
#'                           sigma_subject_intercept = sqrt(0.0333),
#'                           sigma_subject_slope = sqrt(0.0030),
#'                           sigma_error = sqrt(0.0262),
#'                           effect_size = cohend(0.4, standardizer = "slope_SD"))
#'     x <- get_power(p)
#'     round(x$power, 2)
#' }
#'
#' ## Table 1 in Raudenbush & Liu 2001
#' ## NB, it looks like they made an error in column 1.
#' g <- expand.grid(D = 2:8,
#'                  f = c(0.5, 1:6))
#' g$power <- mapply(rauden_liu, D = g$D, f = g$f)
#' tidyr::spread(g, f, power)
#'
#'
#' ## Table 3 Table 1 in Raudenbush & Liu 2001
#' g <- expand.grid(n = seq(100, 800, by = 100),
#'                  D = 4,
#'                  f = c(0.5, 1:6))
#' g$power <- mapply(rauden_liu, n = g$n, f = g$f, D = g$D)
#' tidyr::spread(g, n, power)
#'
cohend <- function(ES, standardizer = "pretest_SD", treatment = "control") {
    if(length(standardizer) != 1) stop("Length of 'standardizer' must be equal to 1", call. = FALSE)
    if(!standardizer %in% c("pretest_SD",
                           "posttest_SD",
                           "slope_SD")) stop("Wrong 'standardizer', allowed options are: 'pretest_SD', 'posttest_SD', and 'slope_SD'.", call. = FALSE)
    if(!treatment %in% c("treatment", "control")) stop("Wrong 'treatment', allowed options are: 'treatment' or 'control'")

    vapply(ES, .cohend,
           standardizer = standardizer,
           treatment = treatment,
           list(1))
}
.cohend <- function(ES, standardizer, treatment) {
    if(length(ES) != 1) stop("Length of ES is not equal to 1.", call. = FALSE)
    # return a ES func
    # that get_slope_diff can evaluate
    if(standardizer == "pretest_SD") {
        f <- calc_slope_from_d(ES, time = "pre", treatment = treatment)
    } else if(standardizer == "posttest_SD") {
        f <- calc_slope_from_d(ES, time = "post", treatment = treatment)
    } else if(standardizer == "slope_SD") {
        f <- function(paras) {
            p <- NA_to_zero(paras)
            p <- prepare_paras(p)
            if(treatment == "control") {
                with(p$control, ES * sqrt(sigma_subject_slope^2 + sigma_cluster_slope^2) * T_end)
            } else {
                with(p$treatment, ES * sqrt(sigma_subject_slope^2 + sigma_cluster_slope^2) * T_end)
            }
        }
    }
    get <- function() {
        list("ES" = ES,
             "standardizer" = standardizer,
             "treatment" = treatment)
    }
    x <- list("set" = f,
              "get" = get)
    class(x) <- append(class(x), "plcp_cohend")

    list(x)

}
# cohend
calc_slope_from_d <- function(ES, time, treatment) {
    function(paras) {
        paras <- NA_to_zero(paras)
        SD <- get_sds(paras, treatment = treatment)
        ind <- ifelse(time == "post", nrow(SD), 1)
        SD <- SD[ind, "SD_with_random_slopes"]
        ES * SD
    }
}
get_effect_size <- function(object) {
        UseMethod("get_effect_size")
}

get_effect_size.plcp <- function(object) {
    ES <- object$effect_size
    if(inherits(ES[[1]], "plcp_cohend")) {
        out <- ES[[1]]$get()
    } else {
        ES <- unlist(ES)
        out <- list("ES" = ES, "standardizer" = "raw", treatment = "")
    }

    out
}
get_effect_size.plcp_multi <- function(object) {
    x <- lapply(1:nrow(object), function(i) {
        data.frame(get_effect_size.plcp(object[i,]))
        }
    )
    x <- do.call(rbind, x)
    x$standardizer <- as.character(x$standardizer)
    x
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


        x$control <- deparse_n2(x$control)
        x$treatment <- deparse_n2(x$treatment)
        tx <- print_per_treatment_(1, x, n2 = TRUE)
        cc <- print_per_treatment_(2, x, n2 = TRUE)
        data.frame(treatment_lab = tx$lab,
                   treatment = tx$num,
                   control_lab = cc$lab,
                   control = cc$num,
                   stringsAsFactors = FALSE)
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


    dropout <- lapply(1:nrow(object), function(i) get_dropout_post(object[i, ]))
    dropout <- do.call(rbind, dropout)
    if(all(dropout$control == dropout$treatment)) {
        out$dropout <- dropout$treatment
    } else {
        out$dropout_tx <- dropout$treatment
        out$dropout_cc <- dropout$control
        out$dropout <- NA
    }

    per_tx_n2 <- vapply(seq_along(object$n2), function(i) is.per_treatment(object$n2[i]), logical(1))
    if(all(n2$treatment_lab == n2$control_lab) & !any(per_tx_n2)) {
        out$n2_lab <- n2$treatment_lab
        out$n2 <- n2$treatment
    } else {
        out$n2_tx_lab <- n2$treatment_lab
        out$n2_tx <- n2$treatment
        out$n2_cc_lab <- n2$control_lab
        out$n2_cc <- n2$control
        #out$n2_tx_lab <- n2$treatment_lab
        #out$n2_cc_lab <- n2$control_lab
        out$n2 <- NULL
    }

    unequal_clust <- lapply(seq_along(object$n2), function(i) is.unequal_clusters(object$n2[i]))
    unequal_clust <- unlist(unequal_clust)

    if(any(unequal_clust)) {
        out$n3 <- NULL
    } else {
        if(all(n3$treatment == n3$control)) {
            out$n3 <- n3$treatment
        } else {
            out$n3_tx <- n3$treatment
            out$n3_cc <- n3$control
            out$n3 <- NULL
        }
    }

    out_dense <- out
    out$icc_pre_cluster <- round(object$icc_pre_cluster, digits)
    out$icc_pre_subject <- round(object$icc_pre_subject, digits)
    out$icc_slope <- round(object$icc_slope, digits)
    out$var_ratio <- round(object$var_ratio, digits)
    ES <- get_effect_size(object)
    out$effect_size <- ES$ES
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
get_multi_title.plcp_mixed <- function(object) {
    "# Multi-study setup (mixed)"
}
get_multi_title.plcp_2lvl <- function(object) {
    "# Multi-study setup (two-level)"
}
get_multi_title.plcp_3lvl <- function(object) {
   "# Multi-study setup (three-level)"
}

select_setup_cols <- function(x) {
    cols <- c("n1",
              "n2_lab", "n2_tx_lab", "n2_cc_lab",
              "dropout", "dropout_tx", "dropout_cc",
              "icc_pre_subject", "icc_pre_cluster", "icc_slope", "var_ratio", "effect_size")
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

    out <- out[, select_setup_cols(out)]
    colnames(out) <- gsub("_lab", "", colnames(out))
    print(out)
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
    attr(n2, "per_treatment") <- FALSE
    n2
}

prepare_paras <- function(paras) {
    paras_tx <- paras
    per_tx_n2 <- FALSE
    if (is.per_treatment(paras$n3)) {
        n3_tx <- paras$n3[[1]]$treatment
        n3_cc <- paras$n3[[1]]$control
        paras$n3 <- n3_cc
        paras_tx$n3 <- n3_tx

       attr(paras$n3, "per_treatment") <- TRUE
       attr(paras_tx$n3, "per_treatment") <- TRUE
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

        per_tx_n2 <- TRUE

    } else {
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

        attr(paras$dropout, "per_treatment") <- TRUE
        attr(paras_tx$dropout, "per_treatment") <- TRUE
    }

    if(length(paras_tx$n2) == 1) {
        paras_tx$n2 <- rep(paras_tx$n2, paras_tx$n3)
    }
    if(length(paras$n2) == 1) {
        paras$n2 <- rep(paras$n2, paras$n3)
    }


    attr(paras$n2, "per_treatment") <- per_tx_n2
    attr(paras_tx$n2, "per_treatment") <- per_tx_n2


    out <- list(control = paras,
         treatment = paras_tx,
         prepared = TRUE)

    class(out) <- class(paras)
    out
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
#' \code{rnorm} or \code{rpois}, then cluster sizes  (and possibly the number of clusters),
#' will be treated as a random variable. The expected power is then reported by averaging over
#' multiple realizations of the random variables.
#'
#' Unless \code{per_treatment} is used, then the same realization of random cluster sizes
#' will be used in both groups. To use independent realizations from the same distribution for
#' each treatment group, simply combine the \code{unequal_clusters} with \code{per_treatment}.
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
#' # Poisson distributed cluster sizes, same in both groups
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
#' # Independent draws from same dist
#' n2 <- unequal_clusters(func = rpois(n = 5, lambda = 5))
#' p <- study_parameters(n1 = 11,
#'                       n2 = per_treatment(n2, n2),
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

    # suport legacy argument 'cohend'
    if("cohend" %in% names(new_args)) {
        new_args$effect_size <- cohend(new_args$cohend)
        new_args$cohend <- NULL
    }
    new <- check_new_argument(args, new_args)
    if(length(new) > 0) stop(paste0("Updating new arguments is not yet implemented. '", new, "' was not used in original call."), call. = FALSE)
    for(i in seq_along(new_args)) {
        args[[names(new_args[i])]] <- new_args[[i]]
    }

    paras <- do.call(study_parameters, args)

    paras
}
#' @method update plcp_multi
#' @export
update.plcp_multi <- function(object, ...) {
    if("plcp_filtered" %in% class(object)) stop("Object is a subset. Update currently only works with the full object.", call. = FALSE)
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
            } else {
                attrib <- attributes(n2_cc)
                n2_cc <- sum(n2_cc)
                attributes(n2_cc) <- attrib
                }
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
