#' Calculate power for two- and three-level models with missing data.
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param df Either "between" or, "satterth" for Satterthwaite's DF approximation.
#' Also accepts a \code{numeric} value which will be used as DF.
#' @param alpha The alpha level, defaults to 0.05.
#' @param progress \code{logical}; displays a progress bar when > 1 power analysis
#' is performed.
#' @param R An \code{integer} indicating how many realizations to base power on.
#' Useful when dropout or cluster sizes are sampled (i.e. are random variables).
#' @param cores An \code{integer} indicating how many CPU cores to use.
#'
#' @param ... Other potential arguments; currently used to pass progress bar from
#'  Shiny
#'
#' @details
#'
#' \bold{Calculation of the standard errors}
#'
#' Designs with equal cluster sizes, and with no missing data, uses standard closed form equations to
#' calculate standard errors. Designs with missing data or unequal cluster sizes uses more
#' computationally intensive linear algebra solutions.
#'
#' To see a more detailed explanation of the calculations, type
#' \code{vignette("technical", package = "powerlmm")}.
#'
#' \bold{Degrees of freedom}
#'
#' Power is calculated using the \emph{t} distribution with non-centrality parameter \eqn{b/se},
#' and \emph{dfs} are either based on a the between-subjects or between-cluster \emph{dfs}, or using Satterthwaite's approximation.
#' For the "between" method, \eqn{N_3 - 2} is used for three-level models, and \eqn{N_2 - 2} for two-level models,
#' where \eqn{N_3} and \eqn{N_2} is the total number of clusters and subjects in both arms.
#'
#' \bold{N.B} Satterthwaite's method will be RAM and CPU intensive for large sample sizes.
#' The computation time will depend mostly on \code{n1} and \code{n2}. For instance, for a fully nested model with
#' \code{n1 = 10}, \code{n2 = 100}, \code{n3 = 4}, computations will likely take 30-60 seconds.
#'
#' \bold{Cluster sizes or dropout pattern that are random (sampled)}
#'
#' If \code{deterministic_dropout = FALSE} the proportion that dropout at each time point will be sampled
#' from a multinomial distribution. However, if it is \code{TRUE}, the proportion of subjects that dropout will be non-random,
#' but which subjects dropout will still be random. Both scenarios often lead to small variations in the estimated power. Moreover,
#' using cluster sizes that are random, \code{unequal_clusters(func = ...)}, can lead to large variations in power
#' for a single realization of cluster sizes. In both scenarios the expected power can be calculated by repeatedly recalculating
#' power for different new realizations of the random variables. This is done be using the argument \code{R} -- power, sample size, and DFs,
#' is then reported by averaging over the \code{R} realizations.
#'
#' If power varies over the \code{R} realization then the Monte Carlo SE is also reported.
#' The SE is based on the normal approximation, i.e. sd(power_i)/sqrt(R).
#'
#' @seealso \code{\link{study_parameters}}, \code{\link{simulate.plcp}}, \code{\link{get_power_table}}
#'
#' @export
#'
#' @return a \code{list} or \code{data.frame} depending if power is calculated for a
#' single set of parameters or a combination of multiple values. Has class
#' \code{plcp_power_3lvl} for fully- and partially nested three-level designs,
#'  and class \code{plcp_power_2lvl} for two-level designs.
#'
#' @examples
#' # Two-level model
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 40,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           var_ratio = 0.02,
#'                           cohend = -0.8)
#'
#' get_power(paras)
#'
#' # With missing data
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 40,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           var_ratio = 0.02,
#'                           dropout = dropout_weibull(0.3, 2),
#'                           cohend = -0.8)
#'
#'
#' get_power(paras)
#'
#'
#' # Three-level model
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 5,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = 0.05,
#'                           var_ratio = 0.02,
#'                           cohend = -0.8)
#'
#' get_power(paras)
#'
#' # With missing data
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 10,
#'                           n3 = 5,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = 0.05,
#'                           var_ratio = 0.02,
#'                           dropout = dropout_weibull(0.3, 2),
#'                           cohend = -0.8)
#'
#' get_power(paras)
#'
#' # Satterthwaite DFs
#' get_power(paras, df = "satterthwaite")
#' @importFrom parallel makeCluster parLapply stopCluster
get_power <- function(object, df = "between", alpha = 0.05, progress = TRUE, R = 1L, cores = 1L, ...) {
    UseMethod("get_power")
}


# print -------------------------------------------------------------------

#' Print method for three-level \code{get_power}
#'
#' @param x An object of class \code{plcp_power_3lvl}.
#' @param ... Optional arguments
#' @method print plcp_power_3lvl
#' @export
print.plcp_power_3lvl <- function(x, ...) {
   .p <- x

   partially_nested <- .p$paras$partially_nested
   MCSE <- 0
   if(.p$R > 0) {
       tot_n <- as.data.frame(.p$tot_n)
       tot_n <- data.frame(control = mean(unlist(tot_n$control)),
                           treatment = mean(unlist(tot_n$treatment)),
                           total = mean(unlist(tot_n$total)))
       n3 <- as.data.frame(do.call(rbind, .p$n3))

       .p$paras$n3 <- per_treatment(mean(unlist(n3$control)),
                                    mean(unlist(n3$treatment)))
       x <- prepare_print_plcp_3lvl(.p$paras)
       width <- attr(x, "width")

      #.p$paras$n2 <- per_treatment(mean(unlist(tot_n$control)),
      #                                mean(unlist(tot_n$treatment)))

       x$total_n <- print_per_treatment(tot_n, width = width)
       if(.p$R > 1) {
           MCSE <- get_monte_carlo_se_gaussian(unlist(.p$power_list))
       }
       .p$power <- mean(unlist(.p$power))
       #x$tot_n <- mean(unlist(.p$tot_n))
       .p$df <- mean(unlist(.p$df))
       #x$se <- mean(unlist(x$se))
   } else x <- prepare_print_plcp_3lvl(.p$paras)


   x$method <- "Power Analyis for Longitudinal Linear Mixed-Effects Models (three-level)\n                  with missing data and unbalanced designs"
   x$df <- .p$df
   x$alpha <- .p$alpha
   if(MCSE > 0) {
       x$power <- paste0(round(unlist(.p$power) * 100, 0), "%", " (MCSE: ", round(MCSE*100), "%)")

   } else {
       x$power <- paste0(round(unlist(.p$power) * 100, 0), "%")
   }

   if(!is.null(x$note) && x$note == "n2 is randomly sampled") {
       txt <- "n2 is randomly sampled"
       x$note <- paste0(txt, ". Values are the mean from R = ", .p$R, " realizations.")
   }
    if(partially_nested) {
        if(is.null(x$note)) {
            x$note <- "Study is partially nested. Clustering only in treatment arm."
        } else {
            x$note <- paste(x$note, "Study is partially nested. Clustering only in treatment arm", sep = "\n      ")
        }

    }
    print(x, ...)

    if(partially_nested & !.p$satterth) {
        message("N.B: Satterthwaite dfs are recommended for partially-nested models, or calculate power with 'simulate.plcp'")
    }
}

#' Print method for two-level \code{get_power}
#'
#' @param x An object of class \code{plcp_power_2lvl}.
#' @param ... Optional arguments
#' @method print plcp_power_2lvl
#' @export
print.plcp_power_2lvl <- function(x, ...) {
    .p <- x
    if(.p$R > 0) {
        tot_n <- .p$tot_n
        tot_n <- data.frame(control = mean(unlist(tot_n$control)),
                            treatment = mean(unlist(tot_n$treatment)),
                            total = mean(unlist(tot_n$total)))

        .p$paras$n2 <- per_treatment(mean(unlist(tot_n$control)),
                                     mean(unlist(tot_n$treatment)))
        x <- prepare_print_plcp_2lvl(.p$paras)
        width <- attr(x, "width")

        #.p$paras$n2 <- per_treatment(mean(unlist(tot_n$control)),
        #                                mean(unlist(tot_n$treatment)))

        #x$total_n <- print_per_treatment(tot_n, width = width)

        .p$power <- mean(unlist(.p$power))
        #x$tot_n <- mean(unlist(.p$tot_n))
        .p$df <- mean(unlist(.p$df))
        #x$se <- mean(unlist(x$se))
    }
    x$df <- .p$df
    x$alpha <- .p$alpha
    x$power <- paste(round(.p$power * 100, 0), "%")
    x$method <- "Power Analysis for Longitudinal Linear Mixed-Effects Models\n            with missing data and unbalanced designs"
    if(.p$R > 1) x$note <- paste0("Sample size is random. Values are the mean from R = ", .p$R, " realizations.")
    print(x)

}




# lmer formual ------------------------------------------------------------


#' Create an lmer formula based on a \code{\link{study_parameters}}-object
#'
#' @param object A \code{\link{study_parameters}}-object containing one study design
#' @param ... Unused, optional arguments.
#' @details
#'
#' The lme4 formula will correspond to the model implied by the specified parameters in
#' the \code{\link{study_parameters}}-object. Thus, if e.g. \code{cor_subject} is \code{NA} or \code{NULL} the
#' corresponding term is removed from the lmer formula. Parameters that are 0 are retained.
#'
#' For crossed design if all correlations involving an effect is NA then that random effect will be modeled as uncorelated.
#'
#'
#' Currently only objects with one study design are supported, i.e. objects with class \code{plcp},
#' and not \code{plcp_multi}; \code{data.frame} with multiple designs are currently not supported.
#'
#' @return A \code{character} vector with lmer formula syntax.
#' @export
create_lmer_formula <- function(object, ...) {
    UseMethod("create_lmer_formula")
}

#' @export
create_lmer_formula.plcp_multi <- function(object, n = 1, ...) {
    if(n > nrow(object)) stop("Row does not exist, 'n' is too large.")
    create_lmer_formula(as.plcp(object[n, ]))
}

#' @export
create_lmer_formula.plcp <- function(object, n = NULL, ...) {
    NextMethod("create_lmer_formula")
}
create_lmer_formula.plcp_nested <- function(object, n = NULL, ...) {
    u0 <- object$sigma_subject_intercept
    u1 <- object$sigma_subject_slope
    u01 <- object$cor_subject

    v0 <- object$sigma_cluster_intercept
    v1 <- object$sigma_cluster_slope
    v01 <- object$cor_cluster

    f0 <- "y ~ time*treatment"
    lvl2 <- make_random_formula(u0, u01, u1, term = "subject")
    if("plcp_2lvl" %in% class(object)) {
        f <- paste(f0, lvl2, sep  = " + ")
    } else if("plcp_3lvl" %in% class(object)) {
        if(object$partially_nested) {
            lvl3 <- make_random_formula_pn(v0, v01, v1)
        } else {
            lvl3 <- make_random_formula(v0, v01, v1, term = "cluster")
        }
        f <-  paste(f0, lvl2, lvl3, sep  = " + ")
    }
    if(object$family == "gaussian") {
        attr(f, "fit_func") <- "lmer"

    } else if(object$family %in% c("binomial")) {
        attr(f, "fit_func") <- "glmer"
        attr(f, "family") <- binomial("logit")
    } else if(object$family %in% c("poisson")) {
        attr(f, "fit_func") <- "glmer"
        attr(f, "family") <- poisson("log")
    } else if(object$family %in% c("gamma")) {
        attr(f, "fit_func") <- "glmer"
        attr(f, "family") <- Gamma("log")
    }

    f
}

get_pars_short_name <- function(object) {
    pars <- list()
    pars["u0"] <- object$sigma_subject_intercept
    pars["u1"] <- object$sigma_subject_slope
    pars["u01"] <- object$cor_subject
    pars["v0"] <- object$sigma_cluster_intercept
    pars["v1"] <- object$sigma_cluster_slope
    pars["v2"] <- object$sigma_cluster_intercept_crossed
    pars["v3"] <- object$sigma_cluster_slope_crossed
    pars["v01"] <- object$cor_cluster_intercept_slope
    pars["v02"] <- object$cor_cluster_intercept_intercept_tx
    pars["v03"] <- object$cor_cluster_intercept_slope_tx
    pars["v12"] <- object$cor_cluster_slope_intercept_tx
    pars["v13"] <- object$cor_cluster_slope_slope_tx
    pars["v23"] <- object$cor_cluster_intercept_tx_slope_tx

    pars
}

create_lmer_formula.plcp_crossed <- function(object, n = NULL, ...) {

    tmp <- get_pars_short_name(object)
    f0 <- "y ~ time*treatment"
    lvl2 <- with(tmp,
                 make_random_formula(u0, u01, u1,
                                     term = "subject")
                 )
    lvl3 <- with(tmp,
                 make_random_formula_crossed(v0, v1, v2, v3,
                                             v01, v02, v03, v12, v13, v23)
                 )
    f <-  paste(f0, lvl2, lvl3, sep  = " + ")

    f
}
make_random_formula <- function(x0, x01, x1, term) {
    if(!is.na(x0) & is.na(x1)) {
        f <- "(1 | g)"
    } else if(is.na(x0) & !is.na(x1)) {
        f <- "(0 + time | g)"
    } else if(!is.na(x0)  & !is.na(x1) & !is.na(x01)) {
        f <- "(1 + time | g)"
    } else if(!is.na(x0) & !is.na(x1) & is.na(x01)) {
        f <- "(1 + time || g)"
    }
    f <- gsub("g", term, f)
    f
}
make_random_formula_pn <- function(x0, x01, x1) {
    if(!is.na(x0) & is.na(x1)) {
        f <- "(0 + treatment | cluster)"
    } else if(is.na(x0) & !is.na(x1)) {
        f <- "(0 + treatment:time | cluster)"
    } else if(!is.na(x0) & !is.na(x1) & !is.na(x01)) {
        f <- "(0 + treatment + treatment:time | cluster)"
    } else if(!is.na(x0) & !is.na(x1) & is.na(x01)) {
        f <- "(0 + treatment + treatment:time || cluster)"
    }
    f
}

## separate crossed slopes that should be correlated or independent
make_re_term <- function(v_i, corr, term) {
    cors <- NULL
    indep <- NULL

    if(!is.na(v_i)) {
        if(corr) cors <- term else indep <- term
    }
    list(cors = cors,
         indep = indep
         )
}
make_cor_formula <- function(terms, subset = "cors") {
    terms <- unlist(terms[,c(subset)])
    if(all(is.na(terms))) return(NULL)
    if(!"1" %in% terms) terms <- c("0", terms)

    slopes <- paste(terms, collapse = " + ")

    if(subset == "cors" | length(terms) == 1) {
        paste("(", slopes, " | cluster)", sep = "")
    } else {
        paste("(", slopes, " || cluster)", sep = "")
    }

}

make_random_formula_crossed <- function(v0, v1, v2, v3, v01, v02, v03, v12, v13, v23) {

    # which cors to include
    # only remove cor when all term involving term is NA
    cor0 <- !all(is.na(c(v01, v02, v03)))
    cor1 <- !all(is.na(c(v01, v12, v13)))
    cor2 <- !all(is.na(c(v02, v12, v23)))
    cor3 <- !all(is.na(c(v03, v13, v23)))

    terms <- mapply(make_re_term,
                    v_i = c(v0, v1, v2, v3),
                    corr = c(cor0, cor1, cor2, cor3),
                    term = c("1", "time", "treatment", "time:treatment"),
                    SIMPLIFY = FALSE)
    terms <- do.call(rbind, terms)

    f <- paste(c(make_cor_formula(terms),
               make_cor_formula(terms, "indep")), collapse = " + ")
    f
}



# new power func ---------------------------------------------------------------

gradient <- function (fun, x, delta = 1e-04, ...)
{
    ## lme4:::grad.ctr3
    nx <- length(x)
    Xadd <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) + diag(delta,
                                                               nx)
    Xsub <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) - diag(delta,
                                                               nx)
    fadd <- apply(Xadd, 1, fun, ...)
    fsub <- apply(Xsub, 1, fun, ...)
    (fadd - fsub)/(2 * delta)
}
make_theta_vec <- function(x0sq, x01, x1sq) {

    # deal with negative paras in numerical derivative
    x0sq <- ifelse(x0sq < 0, abs(x0sq), x0sq)
    #x01 <- ifelse(x01 < 0, abs(x01), x01)
    x1sq <- ifelse(x1sq < 0, abs(x1sq), x1sq)

    if((NA_or_zero(x0sq) | NA_or_zero(x1sq)) & is.na(x01)) {
        x <- c(x0sq, x01, x1sq)
        x <- x[!is.na(x)]
        x <- sqrt(x)
    } else  {
        x <- matrix(c(x0sq, x01, x01, x1sq), ncol = 2)
        if(is.na(x01)) {
           x <- sqrt(x)
           x <- x[c(1,4)]
        # } else {
        #     x <- chol(x)
        #     x <- x[c(1,3,4)]
        } else {
            L <- suppressWarnings(chol(x, pivot = TRUE))
            L <- L[, order(attr(L, "pivot"))]
            x <- L[c(1,3,4)]
        }


    }
    x
}


make_theta_vec_new <- function(x0sq, x01, x1sq) {
    ## not used
    ## meant for crossed-designs

    # deal with negative paras in numerical derivative
    x0sq <- ifelse(x0sq < 0, abs(x0sq), x0sq)
    #x01 <- ifelse(x01 < 0, abs(x01), x01)
    x1sq <- ifelse(x1sq < 0, abs(x1sq), x1sq)
    x <- matrix(c(x0sq, x01,
                  x01, x1sq),
                ncol = 2)
    x0 <- x
    if(all(is.na(x))) return(numeric(0))
    keep <- which(lower.tri(x, diag = TRUE))
    keep <- keep[!is.na(x[keep])]
    del <- which(is.na(c(x0sq,  x1sq)))
    if(length(del) > 0) x <- x[-del, -del, drop = FALSE]

    x[is.na(x)] <- 0

    zeros <- vapply(1:ncol(x), function(i) all(x[i, i] == 0), logical(1))

    x <- x[!zeros, !zeros, drop = FALSE]
    if(ncol(x) == 1) {
        m <- sqrt(x)
    } else {
        x <- nearPD(x, keepDiag = TRUE)$mat
        m <- t(chol(x))
        }

    full <- diag(rep(0, 2))

    m <- m[lower.tri(m, diag = TRUE)]
    keep_zero <- x0[keep] == 0
    full[keep[!keep_zero]] <- m[m != 0]
    full[keep[keep_zero]] <- 0
    full[keep]

}

make_theta <- function(pars) {
    #p <- make_pars(pars)
    p <- as.list(pars)
    sigma <- sqrt(p$sigma)
    #if(old_vec) {
        lvl2 <- make_theta_vec(p$u0, p$u01, p$u1)/sigma
        lvl3 <- make_theta_vec(p$v0, p$v01, p$v1)/sigma
    # } else {
    #     lvl2 <- make_theta_vec(p$u0, p$u01, p$u1)/sigma
    #     lvl3 <- make_theta_vec(p$v0, p$v01, p$v1)/sigma
    # }

    c(lvl2, lvl3)
}
make_theta_crossed <- function(pars) {

    #p <- make_pars(pars)
    p <- as.list(pars)
    sigma <- sqrt(p$sigma)
    lvl2 <- make_theta_vec(p$u0, p$u01, p$u1)/sigma

    x <- with(p, matrix(c(v0, v01, v02,  v03,
                          v01, v1,  v12, v13,
                          v02, v12, v2,   v23,
                          v03, v13, v23,  v3), ncol = 4))

    x0 <- x
    if(all(is.na(x))) return(numeric(0))
    keep <- which(lower.tri(x, diag = TRUE))
    keep <- keep[!is.na(x[keep])]
    del <- with(p,
                which(is.na(c(v0, v1, v2, v3)))
    )
    if(length(del) > 0) x <- x[-del, -del, drop = FALSE]

    x[is.na(x)] <- 0

    zeros <- vapply(1:ncol(x), function(i) all(x[i, i] == 0), logical(1))

    x <- x[!zeros, !zeros, drop = FALSE]
    if(ncol(x) == 1) {
        m <- sqrt(x)
    } else {
        x <- nearPD(x, keepDiag = TRUE)$mat
        m <- t(chol(x))
    }

    full <- diag(rep(0, 2))

    m <- m[lower.tri(m, diag = TRUE)]
    keep_zero <- x0[keep] == 0
    full[keep[!keep_zero]] <- m[m != 0]
    full[keep[keep_zero]] <- 0

    c(lvl2,
      full[keep]/sigma)
}
varb_func <- function(para, X, Zt, L0, Lambdat, Lind, crossed = FALSE) {
    ## adapted from lme4PureR
    ind <- which(!is.na(para))
    pars <- as.list(para)
    function(x = NULL, Lc) {
        if(!is.null(x)) pars[ind] <- x
        if(crossed) {
            theta <- make_theta_crossed(pars)
        } else {
            theta <- make_theta(pars)
        }

        sigma2 <- pars$sigma
        Lambdat@x <- theta[Lind]
        L0 <- Matrix::update(L0, Lambdat %*% Zt, mult = 1)
        XtX <- crossprod(X)
        ZtX <- Zt %*% X
        RZX <- Matrix::solve(L0, Matrix::solve(L0, Lambdat %*% ZtX, system = "P"),
                             system = "L")
        RXtRX <- as(XtX - crossprod(RZX), "dpoMatrix")

        t(Lc) %*% (sigma2 * solve(RXtRX)) %*% Lc
    }
}
setup_power_calc <- function(object, d, f) {
    UseMethod("setup_power_calc")
}


setup_power_calc.plcp_nested <- function(object, d, f) {
    u0 <- object$sigma_subject_intercept
    u1 <- object$sigma_subject_slope
    cor_subject <- object$cor_subject
    u01 <- u0 * u1 * cor_subject
    v0 <- object$sigma_cluster_intercept
    v1 <- object$sigma_cluster_slope
    v01 <- v0 * v1 * object$cor_cluster
    sigma <- object$sigma_error
    sigma2 <- sigma^2
    pars <- c("u0" = u0^2, "u01" = u01, "u1" = u1^2,
              "v0" = v0^2, "v01" = v01, "v1" = v1^2,
              "sigma" = sigma^2)
    theta <- make_theta(pars)

    X <- f$X
    Lambdat <- f$reTrms$Lambdat
    Lind <- f$reTrms$Lind
    Zt <- f$reTrms$Zt
    Lambdat@x <- theta[Lind]
    L0 <- Matrix::Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult = 1)

    list("pars" = pars,
         "theta" = theta,
         "X" = X,
         "Zt" = Zt,
         "Lambdat" = Lambdat,
         "L0" = L0,
         "Lind" = Lind)

}
setup_power_calc.plcp_crossed <- function(object, d, f) {

    pars <- get_pars_short_name(object)
    pars <- with(pars,
                 list(u0 = u0^2,
                 u1 = u1^2,
                 u01 = u0 * u1 * u01,
                 v0 = v0^2,
                 v1 = v1^2,
                 v2 = v2^2,
                 v3 = v3^2,
                 v01 = v0 * v1 * v01,
                 v02 = v0 * v2 * v02,
                 v03 = v0 * v3 * v03,
                 v12 = v1 * v2 * v12,
                 v13 = v1 * v3 * v13,
                 v23 = v2 * v3 * v23,
                 sigma = object$sigma_error^2)
                 )

    theta <- make_theta_crossed(pars)

    X <- f$X
    Lambdat <- f$reTrms$Lambdat
    Lind <- f$reTrms$Lind
    Zt <- f$reTrms$Zt
    Lambdat@x <- theta[Lind]
    L0 <- Matrix::Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult = 1)

    list("pars" = pars,
         "theta" = theta,
         "X" = X,
         "Zt" = Zt,
         "Lambdat" = Lambdat,
         "L0" = L0,
         "Lind" = Lind)

}

power_worker <- function(object, df, alpha, use_satterth, ...) {
    dots <- list(...)
    use_matrix_manual <- ifelse(is.null(dots$use_matrix), FALSE, dots$use_matrix)
    crossed <- inherits(object, "plcp_crossed")
    function(i = NULL) {
        use_matrix_se <- is.unequal_clusters(object$n2) | is.list(object$dropout) | use_satterth
        prepped <- prepare_paras(object)
        if(use_matrix_se || use_matrix_manual) {
            print("useMatrix")
            d <- simulate_data(prepped)
            f <- lme4::lFormula(formula = create_lmer_formula(object),
                                data = d)

            pc <- setup_power_calc(object = object,
                                   d = d,
                                   f = f)
            pars <- pc$pars
            X <- pc$X
            Zt <- pc$Zt
            L0 <- pc$L0
            Lambdat <- pc$Lambdat
            Lind <- pc$Lind

            varb <- varb_func(para = pars,
                              X = X,
                              Zt = Zt,
                              L0 = L0,
                              Lambdat = Lambdat,
                              Lind = Lind,
                              crossed = crossed)
            Phi <- varb(Lc = diag(4))
            se <- sqrt(Phi[4, 4])
            calc_type <- "matrix"
        } else {
            se <- get_se_classic(prepped)
            calc_type <- "classic"
        }

        if(use_satterth) {
            df <- get_satterth_df(prepped,
                                  d = d,
                                  pars = pars,
                                  Lambdat = Lambdat,
                                  X = X,
                                  Zt = Zt,
                                  L0 = L0,
                                  Phi = Phi,
                                  varb = varb)
        } else if(df == "between") {
            df <- get_balanced_df(object)
        } else if(is.numeric(df)) df <- df

        # power
        slope_diff <- get_slope_diff(object)/object$T_end
        lambda <- slope_diff / se

        power <- pt(qt(1-alpha/2, df = df), df = df, ncp = lambda, lower.tail = FALSE) +
            pt(qt(alpha/2, df = df), df = df, ncp = lambda)
        tot_n <- list(control = sum(unlist(prepped$control$n2)),
                      treatment = sum(unlist(prepped$treatment$n2)))
        n2 <- list(control = prepped$control$n2,
                   treatment = prepped$treatment$n2)
        n3 <- list(control = prepped$control$n3,
                   treatment = prepped$treatment$n3)

        tot_n$total <- tot_n$control + tot_n$treatment
        list(power = power,
             df = df,
             n2 = n2,
             n3 = n3,
             tot_n = tot_n,
             se = se,
             calc_type = calc_type)
    }
}


unnest_tot_n <- function(x) {
    x <- do.call(rbind, x)
    x <- as.data.frame(x)
    for(i in 1:ncol(x)) {
        x[ ,i] <- unlist(x[,i])
    }
    x
}
unnest_n2 <- function(x) {
    x <- do.call(rbind, x)
    x <- as.data.frame(x)
    tmp <- vector("list", ncol(x))
    for(i in 1:ncol(x)) {
        tmp[[i]] <-  do.call(rbind, x[,i])
    }
    names(tmp) <- colnames(x)
    tmp
}

#' @export
get_power.plcp <- function(object, df = "between", alpha = 0.05, progress = TRUE, R = 1L, cores = 1L, ...) {
    if(R == 1) progress <- FALSE
    dots <- list(...)
    cl <- dots$cl
   # if(is.null(d)) d <- simulate_data(object)
    use_satterth <- (df == "satterthwaite" | df == "satterth")

    power_fun <- power_worker(object = object,
                              df = df,
                              alpha = alpha,
                              use_satterth = use_satterth,
                              ...)
    if(cores > 1) {
        if(is.null(cl)) {
            cl <- makeCluster(getOption("cl.cores", min(R, cores)))
            stop_cluster <- TRUE
        } else stop_cluster <- FALSE
        parallel::clusterEvalQ(cl, expr =
                                   suppressPackageStartupMessages(require(powerlmm, quietly = TRUE)))
        #clusterExport(cl = cl, envir = environment())
        tmp <- parLapply(cl, X = 1:R, power_fun)

        if(stop_cluster) stopCluster(cl)
    } else {
        if(progress) pb <- txtProgressBar(style = 3, min = 1, max = R)
        tmp <- vector(mode = "list", length = R)
        for(i in 1:R) {
            if(progress) setTxtProgressBar(pb, i)
            tmp[[i]] <- power_fun(i)
        }
        if(progress) close(pb)
    }
    tmp <- as.data.frame(do.call(rbind, tmp))

    power <- unlist(tmp$power)
    power_list <- power
    power <- mean(power)
    df <- unlist(tmp$df)
    df_list <- df
    df <- mean(df)
    se <- unlist(tmp$se)
    se_list <- se
    se <- mean(se)
    calc_type <- tmp$calc_type
    n2 <- unnest_n2(tmp$n2)
    tot_n <- unnest_tot_n(tmp$tot_n)


    n3 <- tmp$n3

    out <- list(power = power,
                power_list = power_list,
                df = df,
                df_list = df_list,
                satterth = use_satterth,
                se = se,
                se_list = se_list,
                paras = object,
                alpha = alpha,
                calc_type = calc_type,
                tot_n = tot_n,
                n2 = n2,
                n3 = n3,
                R = R)

    if("plcp_2lvl" %in% class(object))  class(out) <- append(class(out), "plcp_power_2lvl")
    if("plcp_3lvl" %in% class(object))  class(out) <- append(class(out), "plcp_power_3lvl")

    out
}


multi_power_worker <- function(object, df, alpha, R, ...) {

    out <- get_power.plcp(object,
                          df = df,
                          alpha = alpha,
                          R = R,
                          ...)
    tot_n <- out$tot_n
    tot_n <- as.data.frame(tot_n)
    power <- unlist(out$power)
    power_list <- unlist(out$power_list)
    se <- unlist(out$se)
    se_list <- out$se_list
    df <- unlist(out$df)
    df_list <- out$df_list
    out <- list(power = power,
                power_SD = sd(power_list),
                power_list = power_list,
                df = df,
                df_list = df_list,
                tot_n = tot_n,
                se = se,
                se_list = se_list,
                n2_list = out$n2)
}
loop_power <- function(object, df, alpha, nr, progress, progress_inner = FALSE, R, ...) {
    if(progress) pb <- txtProgressBar(style = 3, min = 0, max = nr)
    x <- vector(mode = "list", length = nr)
    for(i in 1:nrow(object)) {
        p <- as.plcp(object[i,])
        out <- multi_power_worker(object = p,
                                df = df,
                                alpha = alpha,
                                R = R,
                                progress = progress_inner,
                                ...)
        if(progress) setTxtProgressBar(pb, i)
        x[[i]] <- out
    }
    if(progress) close(pb)

    x
}

#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar

get_power.plcp_multi <- function(object, df = "between", alpha = 0.05, progress = TRUE, R = 1, cores = 1, ...) {
    dots <- list(...)
    if (is.function(dots$updateProgress)) {
        dots$updateProgress()
    }
    nr <- nrow(object)
    if(cores == 1) {
        progress_inner <- progress
        progress <- FALSE
        x <- loop_power(object = object,
                        df = df,
                        alpha = alpha,
                        nr = nr,
                        progress = progress,
                        progress_inner = progress_inner,
                        R = R)

    } else {
        if(is.null(dots$cl)) {
            cl <- parallel::makeCluster(min(cores, max(nr, R)))
            on.exit(parallel::stopCluster(cl))
        } else cl <- dots$cl
        parallel::clusterEvalQ(cl, expr =
                                   suppressPackageStartupMessages(require(powerlmm, quietly = TRUE)))
        parallel::clusterExport(cl, "multi_power_worker", envir = parent.env(environment()))
        #parallel::clusterExport(cl, "R", envir = environment())

        if(R > nr) {
            x <- loop_power(object = object,
                            df = df,
                            alpha = alpha,
                            nr = nr,
                            progress = progress,
                            cl = cl,
                            R = R,
                            cores = cores)
        } else {
            p <- vector("list", nr)
            for(i in 1:nrow(object)) {
                p[[i]] <- as.plcp(object[i,])
            }
            x <- parLapply(cl,
                           X = p,
                           fun = multi_power_worker,
                           df = df,
                           alpha = alpha,
                           R = R)
        }
    }

    x <- do.call(rbind, x)
    x <- x[, c("power",
               "power_SD",
               "tot_n",
               "power_list",
               "df",
               "se",
               "n2_list"), drop = FALSE]

    out_dense <- prepare_multi_power_out(object,
                                         x = x,
                                         R = R,
                                         alpha = alpha,
                                         df = df)

    out_dense
}

prepare_multi_power_out <- function(object, x, R, alpha, df) {
    .x <- x
    x <- cbind(object, x)
    prep <- prepare_multi_setup(object)
    out <- prep$out
    per_treatment <- all(colnames(out) != "n2")
    if(per_treatment) {
        out$n2_tx <- truncate_n2(out$n2_tx)
        out$n2_cc <- truncate_n2(out$n2_cc)
    } else {
        out$n2 <- truncate_n2(out$n2)
    }
    per_treatment_n3 <- all(colnames(out) != "n3")
    if(per_treatment_n3) {
        out$n3_tx <- out$n3_tx
        out$n3_cc <- out$n3_cc
    } else {
        out$n3 <- out$n3
    }

    out_dense <- prep$out_dense
    out <- out[, select_setup_cols(out)]
    out$df <- round(unlist(x$df), 2)
    out$power <- paste(round(unlist(x$power) * 100, 1), "%")
    if(R > 1) {
        MCSE <- lapply(x$power_list, get_monte_carlo_se_gaussian)
        out$power_MCSE <- paste(round(unlist(MCSE) * 100, 1), "%")
    }
    out_dense$df <- unlist(x$df)
    out_dense$power <- unlist(x$power)
    out_dense$power_SD <- unlist(x$power_SD)
    out_dense$power_list <- x$power_list
    out_dense$tot_n <- x$tot_n
    out_dense$se <- unlist(x$se)
    #out_dense <- cbind(out_dense, ES)
    out_dense$n2_list <- x$n2_list

    class(out_dense) <- append("plcp_multi_power", class(out_dense))
    attr(out_dense, "out") <- out
    attr(out_dense, "x") <- .x
    attr(out_dense, "object") <- object
    attr(out_dense, "R") <- R
    attr(out_dense, "alpha") <- alpha
    attr(out_dense, "df") <- df

    out_dense
}

#' Print method for \code{get_power}-multi
#'
#' @param x An object of class \code{plcp_multi_power}.
#' @param ... Optional arguments
#' @method print plcp_multi_power
#' @export
print.plcp_multi_power <- function(x, ...) {

     out <- attr(x, "out")
    alpha <- attr(x, "alpha")
    df <- attr(x, "df")
    R <- attr(x, "R")
    out <- as.data.frame(out)
    #cat(get_multi_title(x$object), "\n")
    colnames(out) <- gsub("_lab", "", colnames(out))
    cat("# Power Analysis for Longitudinal Linear Mixed-Effect Models\n\n")
    print(out, digits = 2)
    cat(paste0("---\n# alpha = ", alpha, "; DFs = ", df, "; R = ", R), "\n")
    invisible(x)
}

#' Subset function for \code{plcp_multi_power}-objects
#'
#' Custom subset function for \code{plcp_multi_power}-object to make it compatible
#' with its print method.
#' @param x A \code{plcp_multi_power}-object.
#' @param i Indicates which rows to subset.
#' @param ... Ignored.
#' @method [ plcp_multi_power
#' @export
`[.plcp_multi_power` <- function(x, i, ...) {

    if(length(i) == 1 && i > nrow(x)) stop("Row number does not exist.", call. = FALSE)
    if(all(!i)) stop("Nothing to print, subset empty.", call. = FALSE)

    x_new <- attr(x, "x")[i, ]
    object <- attr(x, "object")[i, ]

    out_dense <- prepare_multi_power_out(object,
                                         x = x_new,
                                         R = attr(x, "R"),
                                         alpha = attr(x, "alpha"),
                                         df = attr(x, "df"))
    class(out_dense) <-  append(class(out_dense), "plcp_filtered")

    out_dense
}
