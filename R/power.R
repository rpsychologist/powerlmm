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
#' and DFs are either based on a the between-subjects or between-cluster \emph{dfs}, or using Satterthwaite's approximation.
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
            x$note <- "Study is partially-nested. Clustering only in treatment arm. Cohen's d is standardized using the control group's pretest SD."
        } else {
            x$note <- paste(x$note, "Study is partially-nested. Clustering only in treatment arm", sep = "\n      ")
        }

    }
    print(x, ...)

    if(partially_nested & !.p$satterth) {
        message("N.B: Satterthwaite DFs are recommended for partially-nested models, or calculate power with 'simulate.plcp'")
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

    # deal with negative paras in numeric derivative
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
make_theta <- function(pars) {
    #p <- make_pars(pars)
    p <- as.list(pars)
    sigma <- sqrt(p$sigma)
    lvl2 <- make_theta_vec(p$u0, p$u01, p$u1)/sigma
    lvl3 <- make_theta_vec(p$v0, p$v01, p$v1)/sigma
    c(lvl2, lvl3)
}
varb_func <- function(para, X, Zt, L0, Lambdat, Lind) {
    ## adapted from lme4PureR
    ind <- which(!is.na(para))
    pars <- as.list(para)
    function(x = NULL, Lc) {
        if(!is.null(x)) pars[ind] <- x
        theta <- make_theta(pars)
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

setup_power_calc <- function(d, f, object) {
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

power_worker <- function(object, df, alpha, use_satterth) {

    function(i = NULL) {
        use_matrix_se <- is.unequal_clusters(object$n2) | is.list(object$dropout) | use_satterth
        prepped <- prepare_paras(object)
        if(use_matrix_se) {
            d <- simulate_data(prepped)
            f <- lme4::lFormula(formula = create_lmer_formula(object),
                                data = d)

            pc <- setup_power_calc(d, f, object)
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
                              Lind = Lind)
            Phi <- varb(Lc = diag(4))
            se <- sqrt(Phi[4,4])
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
                              use_satterth = use_satterth)
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



# OLD ---------------------------------------------------------------------


get_power_3lvl_old.list <- function(object, ...) {
    paras <- object
    dots <- list(...)
    n1 <- paras$n1

    if(!is.unequal_clusters(paras$n2) & !is.per_treatment(paras$n2)) {
        paras$n2 <- per_treatment(unlist(paras$n2),
                                  unlist(paras$n2))
    }
    if(!is.per_treatment(paras$n3)) {
        paras$n3 <- per_treatment(unlist(paras$n3),
                                  unlist(paras$n3))
    }
    n2 <- paras$n2
    n3 <- paras$n3
    T_end <- paras$T_end

    error <- paras$sigma_error
    u1 <- paras$sigma_subject_slope
    v1 <- paras$sigma_cluster_slope

    slope_diff <- get_slope_diff(paras)/T_end

    res <- get_power_3lvl.paras(n1 = n1,
                                n2 = n2,
                                n3 = n3,
                                T_end = T_end,
                                error = error,
                                u1 = u1,
                                v1 = v1,
                                slope_diff = slope_diff,
                                partially_nested = paras$partially_nested,
                                allocation_ratio = paras$allocation_ratio,
                                dropout = paras$dropout,
                                paras = paras,
                                ...)


    # show progress in shiny
    if (is.function(dots$updateProgress)) {
        dots$updateProgress()
    }

    res <- list(
        power = res$power,
        dropout_cc = list(as.character(res$dropout_cc)),
        se = res$se,
        df = res$df,
        dropout_tx = list(as.character(res$dropout_tx)),
        n1 = paras$n1,
        n2_cc = res$n2_cc,
        n2_tx = res$n2_tx,
        n3_cc = res$n3_cc,
        n3_tx = res$n3_tx,
        allocation_ratio = paras$allocation_ratio,
        tot_n = res$tot_n,
        var_ratio = get_var_ratio(paras),
        icc_slope = get_ICC_slope(paras),
        icc_pre_subjects = get_ICC_pre_subjects(paras),
        icc_pre_clusters = get_ICC_pre_clusters(paras),
        cohend = paras$cohend,
        T_end = paras$T_end,
        partially_nested = res$partially_nested,
        paras = paras)

    class(res) <- append("plcp_power_3lvl", class(res))

    res

}
get_power_2lvl.list <- function(object, ...) {
    paras <- object
    dots <- list(...)
    n1 <- paras$n1

    tmp <- get_tot_n(paras)
    paras$n2 <- per_treatment(tmp$control, tmp$treatment)
    paras$n3 <- per_treatment(1, 1)
    n2 <- paras$n2

    n3 <- paras$n3
    T_end <- paras$T_end

    error <- paras$sigma_error
    u1 <- paras$sigma_subject_slope
    v1 <- 0

    slope_diff <- get_slope_diff(paras)/T_end

    res <- get_power_3lvl.paras(n1 = n1,
                                n2 = n2,
                                n3 = n3,
                                T_end = T_end,
                                error = error,
                                u1 = u1,
                                v1 = v1,
                                slope_diff = slope_diff,
                                partially_nested = paras$partially_nested,
                                allocation_ratio = paras$allocation_ratio,
                                dropout = paras$dropout,
                                paras = paras,
                                ...)


    # show progress in shiny
    if (is.function(dots$updateProgress)) {
        dots$updateProgress()
    }

    res <- list(power = res$power,
                dropout_cc = list(as.character(res$dropout_cc)),
                se = res$se,
                df = res$df,
                dropout_tx = list(as.character(res$dropout_tx)),
                n1 = paras$n1,
                n2_cc = sum(unlist(res$n2_cc)),
                n2_tx = sum(unlist(res$n2_tx)),
                tot_n = res$tot_n,
                allocation_ratio = paras$allocation_ratio,
                var_ratio = get_var_ratio(paras),
                icc_slope = get_ICC_slope(paras),
                icc_pre_subjects = get_ICC_pre_subjects(paras),
                cohend = paras$cohend,
                T_end = paras$T_end,
                paras = paras)

    class(res) <- append("plcp_power_2lvl", class(res))

    res

}

get_power_3lvl.paras <- function(n1,
                                 n2,
                                 n3,
                                 T_end,
                                 error,
                                 u1,
                                 v1,
                                 slope_diff,
                                 partially_nested,
                                 allocation_ratio = 1,
                                 dropout = NULL,
                                 ...) {
    dots <- list(...)

    sx <- Vectorize(var_T)(n1, T_end)


    if(is.unequal_clusters(n2) | is.list(dropout)) {

        res <- get_se_3lvl_matrix(dots$paras)
        se <- res$se
        var_cc <- res$var_cc
        var_tx <- res$var_tx

        dropout_cc <- format_dropout(get_dropout(dots$paras)$control)
        dropout_tx <- format_dropout(get_dropout(dots$paras)$treatment)
    } else {
        n2_tx <- n2[[1]]$treatment
        n2_cc <- get_n2(dots$paras)$control

        n3_tx <- n3[[1]]$treatment
        n3_cc <- n3[[1]]$control
        if(partially_nested) {
            se <- sqrt( (error^2 + n1*u1^2*sx)/(n1*n2_cc*sx) + (error^2 + n1*u1^2*sx + n1*n2_tx*v1^2*sx) / (n1*n2_tx*n3_tx*sx) )
        } else {
            var_cc <- (error^2 + n1*u1^2*sx + n1*n2_cc*v1^2*sx) / (n1*n2_cc*n3_cc*sx)
            var_tx <- (error^2 + n1*u1^2*sx + n1*n2_tx*v1^2*sx) / (n1*n2_tx*n3_tx*sx)
            se <- sqrt(var_cc + var_tx)
        }
        dropout_cc <- 0
        dropout_tx <- 0

    }
    n2_tx <- get_tot_n(dots$paras)$treatment
    n2_cc <- get_tot_n(dots$paras)$control

    n3_cc <- get_n3(dots$paras)$control
    n3_tx <- get_n3(dots$paras)$treatment


    lambda <- slope_diff / se
    if(v1 == 0) {
        df <-  (n2_tx + n2_cc) - 2
    } else if(!partially_nested) {
        df <- (n3_cc + n3_tx) - 2
    } else {
        df <-(2*n3_tx) - 2

    }
    power <- pt(qt(1-0.05/2, df = df), df = df, ncp = lambda, lower.tail = FALSE) +
        pt(qt(0.05/2, df = df), df = df, ncp = lambda)

    # show progress in shiny
    if (is.function(dots$updateProgress)) {
        dots$updateProgress()
    }

    res <- data.frame(n1 = n1,
                      n3_cc = n3_cc,
                      n3_tx = n3_tx,
                      ratio = get_var_ratio(u1=u1, v1=v1, error=error),
                      icc_slope = get_ICC_slope(u1=u1, v1=v1),
                      tot_n = get_tot_n(dots$paras)$total,
                      se = se,
                      df = df,
                      power = power,
                      dropout_tx = dropout_tx,
                      dropout_cc = dropout_cc,
                      partially_nested = partially_nested)

    res$n2_tx <- list(n2_tx)
    res$n2_cc <- list(n2_cc)



    res
}

get_se_classic <- function(object) {

    if(is.null(object$prepared)) {
        p_tx <- prepare_paras(object)
    } else {
        p_tx <- object
        object <- p_tx$treatment
    }

    p_cc <- p_tx$control
    p_tx <- p_tx$treatment

    n1 <- object$n1
    T_end <- object$T_end
    sx <- Vectorize(var_T)(n1, T_end)



    n2_tx <- unique(unlist(p_tx$n2))
    n3_tx <- unlist(p_tx$n3)
    n2_cc <- unique(unlist(p_cc$n2))
    n3_cc <- unlist(p_cc$n3)


    error <- object$sigma_error
    u1 <- object$sigma_subject_slope
    u1[is.na(u1)] <- 0
    v1 <- object$sigma_cluster_slope
    v1[is.na(v1)] <- 0

    if(object$partially_nested) {
        se <- sqrt( (error^2 + n1*u1^2*sx)/(n1*n2_cc*n3_cc*sx) + (error^2 + n1*u1^2*sx + n1*n2_tx*v1^2*sx) / (n1*n2_tx*n3_tx*sx) )
    } else {
        var_cc <- (error^2 + n1*u1^2*sx + n1*n2_cc*v1^2*sx) / (n1*n2_cc*n3_cc*sx)
        var_tx <- (error^2 + n1*u1^2*sx + n1*n2_tx*v1^2*sx) / (n1*n2_tx*n3_tx*sx)
        se <- sqrt(var_cc + var_tx)
    }


    se

}



# Matrix power ------------------------------------------------------------

## Unbalanced
create_Z_block <- function(n2) {
    B <- matrix(c(1, 0, 0, 1), nrow=2, ncol=2)
    A <- array(1, dim = c(n2, 1))
    kronecker(A, B)
}


#' @import Matrix
get_vcov <- function(paras) {
    n1 <- paras$n1
    n2 <- unlist(paras$n2)
    n3 <- paras$n3
    if(length(n2) == 1) {
        n2 <- rep(n2, n3)
    }

    # paras
    u0 <- paras$sigma_subject_intercept
    u1 <- paras$sigma_subject_slope
    u01 <- u0 * u1 * paras$cor_subject
    v0 <- paras$sigma_cluster_intercept
    v1 <- paras$sigma_cluster_slope
    v01 <- v0 * v1 * paras$cor_cluster
    sigma <- paras$sigma_error
    lvl2_re <- matrix(c(u0^2, u01,
                        u01, u1^2), nrow = 2, ncol = 2)
    lvl3_re <- matrix(c(v0^2, v01,
                        v01, v1^2), nrow = 2, ncol = 2)


    ##
    tot_n <- get_tot_n(paras)$control # tx == cc
    A <- Diagonal(tot_n)
    B <- Matrix(c(rep(1, n1), get_time_vector(paras)), ncol = 2, nrow = n1)
    X <- kronecker(A, B)

    # missing
    if(is.list(paras$dropout) | is.function(paras$dropout[[1]])) {
        miss <- dropout_process(unlist(paras$dropout), paras)
        cluster <- create_cluster_index(n2, n3)
        cluster <- rep(cluster, each = n1)
        miss$cluster <- cluster
        miss <- miss[order(miss$id), ]
        X <- X[miss$missing == 0, ]
    }
    Xt <- Matrix::t(X)
    ## random
    I1 <- Diagonal(tot_n)
    lvl2 <- kronecker(I1, lvl2_re)

    I3 <- Diagonal(n3)
    Z <- bdiag(lapply(n2, create_Z_block))

    lvl3 <- Z %*% kronecker(I3, lvl3_re)
    lvl3 <- Matrix::tcrossprod(lvl3, Z)


    # missing data
    if(is.list(paras$dropout) | is.function(paras$dropout[[1]])) {
        ids <- miss[miss$missing == 0, ]

        ids <- lapply(unique(ids$id), function(id) {
            n <- length(ids[ids$id == id, 1])
            data.frame(id = id,
                       n = n)
        })
        ids <- do.call(rbind, ids)
        ids <- ids[ids$n == 1, "id"]

        s_id <- c(2*ids - 1, 2*ids)
        tmp <- Xt %*% X
        if(length(s_id) == 0) {
            tmp <- solve(tmp)
        } else {
            tmp[-s_id, -s_id] <- solve(tmp[-s_id, -s_id])
        }
    } else {
        ids <- NULL
        s_id <- NULL
        tmp <- solve(Xt %*% X)
    }
    XtX_inv <- tmp
    if(length(ids) > 0) {
        lvl2[ids*2, ] <- 0
        lvl2[, ids*2] <- 0

        lvl3[ids*2, ] <- 0
        lvl3[, ids*2] <- 0
    }
    I <- Diagonal(nrow(X))
    A <- sigma^-2 * (I - X %*% XtX_inv %*% Xt)
    B <- sigma^2 * XtX_inv + (lvl2 + lvl3)

    B_inv <- B
    rm(B)

    # deal with subjects with only 1 observations
    if(length(s_id) == 0) {
        B_inv <- solve(B_inv)
    } else {
        tmp <- Matrix::solve(B_inv[-s_id, -s_id])
        tmp <- as.matrix(tmp)
        B_inv <- as.matrix(B_inv)
        B_inv[-s_id, -s_id] <- tmp
        B_inv <- Matrix(B_inv)
        rm(tmp)
        if(length(ids) == 1) {
            B_inv[2*ids-1, 2*ids-1] <-  1/B_inv[2*ids-1, 2*ids-1]
        } else {
            Matrix::diag(B_inv[2*ids-1, 2*ids-1]) <-
                1/Matrix::diag(B_inv[2*ids-1, 2*ids-1])
        }
    }

    C <- X %*% XtX_inv %*% B_inv %*% XtX_inv %*% Xt
    V_inv <- A+C

    W <- create_Z_block(n3)
    var_betas <- solve(Matrix::t(W) %*% Matrix::t(Z) %*% Xt %*%
                           V_inv %*% X %*% Z %*% W)

    var_betas
}


get_se_3lvl_matrix <- function(paras, ...) {
    dots <- list(...)

    unequal_allocation <- is.per_treatment(paras$n2) | is.per_treatment(paras$n3)
    dropout_per_treatment <- is.per_treatment(paras$dropout)

    tmp <- prepare_paras(paras)
    paras <- tmp$control
    paras_tx <- tmp$treatment

    var_grp1 <- get_vcov(paras)

    if(unequal_allocation | dropout_per_treatment | paras$partially_nested) {
        var_grp2 <- get_vcov(paras_tx)
    } else {
        var_grp2 <- var_grp1
    }

    se <- sqrt(var_grp1[2,2] + var_grp2[2,2])

    list(se = se,
         var_cc = var_grp1[2,2],
         var_tx = var_grp2[2,2])
}


