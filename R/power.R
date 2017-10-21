#' Calculate power for two- and three-level models with missing data.
#'
#' @param object An object created by \code{\link{study_parameters}}
#'
#' @param ... Other potential arguments; currently used to pass progress bar from
#'  Shiny
#'
#' @details
#'
#' \bold{Calculations of standard errors}
#' Designs with equal cluster sizes, and with no missing data, uses standard closed form equations to
#' calculate standard errors. Designs with missing data or unequal cluster sizes uses more
#' computationally intensive linear algebra solutions.
#'
#' To see a more detailed explanation of the calculations, type
#' \code{vignette("technical", package = "powerlmm")}.
#'
#' \bold{Degrees of freedom}
#' Power is calculated using the \emph{t} distribution with non-centrality parameter \eqn{d/se}.
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
get_power <- function(object, ...) {
    UseMethod("get_power")
}
## 2 lvl power
get_power_2lvl <- function(object, ...) {
    UseMethod("get_power_2lvl")

}
#' @export
get_power.plcp_2lvl <- function(object, ...) {
    get_power_2lvl(object, ...)
}

get_power_2lvl.data.frame <- function(object, ...) {

    res <- lapply(1:nrow(object), function(i) get_power_2lvl.list(as.plcp(object[i,]), ...))
    res <- do.call(rbind, res)
    res <- as.data.frame(res)
    res

}

#' @export
get_power.plcp_3lvl <- function(object, ...) {
    get_power_3lvl(object, ...)
}

get_power_3lvl <- function(object, ...) {
    UseMethod("get_power_3lvl")

}
get_power_3lvl.data.frame <- function(object, ...) {
    res <- lapply(1:nrow(object), function(i) get_power_3lvl.list(as.plcp(object[i,]), ...))
    res <- do.call(rbind, res)
    res <- as.data.frame(res)
    res
}

get_power_3lvl.list <- function(object, ...) {
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

# print -------------------------------------------------------------------

#' Print method for three-level \code{get_power}
#'
#' @param x An object of class \code{plcp_power_3lvl}.
#' @param ... Optional arguments
#' @method print plcp_power_3lvl
#' @export
print.plcp_power_3lvl <- function(x, ...) {
   .p <- x
   x <- prepare_print_plcp_3lvl(.p$paras)
   x$method <- "Power calculation for longitudinal linear mixed model (three-level)\n                           with missing data and unbalanced designs"
   x$power <- .p$power
    if(.p$partially_nested) {
        x$note <- "Study is partially-nested. Clustering only in treatment arm"
    }

    print(x, digits = 2, ...)

    if(.p$partially_nested) {
        message("N.B: The degrees of freedom for partially nested designs are experimental, see '?get_power'")
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
    x <- prepare_print_plcp_2lvl(.p$paras)
    x$power <- .p$power
    x$method <- "Power calculation for longitudinal linear mixed model\n            with missing data and unbalanced designs"
    print(x, digits = 2)

}


#

# lmer formual ------------------------------------------------------------


#' Create an lmer formula based on a \code{\link{study_parameters}}-object
#'
#' @param object A \code{\link{study_parameters}}-object containing one study design
#' @details
#'
#' The lmer formula will correspond to the model implied by the non-zero parameters in
#' the \code{\link{study_parameters}}-object. Thus, if e.g. \code{cor_subject} is 0 the
#' corresponding term is removed from the lmer formula.
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
create_lmer_formula.plcp_multi <- function(object, n = 1) {
    if(n > nrow(object)) stop("Row does not exist, 'n' is too large.")
    create_lmer_formula(as.plcp(p1[n,]))
}

#' @export
create_lmer_formula.plcp <- function(object, n = NULL) {
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
    if(x0 != 0 & x1 == 0) {
        f <- "(1 | g)"
    } else if(x0 == 0 & x1 != 0) {
        f <- "(0 + time | g)"
    } else if(x0 != 0 & x1 != 0 & x01 != 0) {
        f <- "(1 + time | g)"
    } else if(x0 != 0 & x1 != 0 & x01 == 0) {
        f <- "(1 + time || g)"
    }
    f <- gsub("g", term, f)
    f
}
make_random_formula_pn <- function(x0, x01, x1) {
    if(x0 != 0 & x1 == 0) {
        f <- "(0 + treatment | cluster)"
    } else if(x0 == 0 & x1 != 0) {
        f <- "(0 + treatment:time | cluster)"
    } else if(x0 != 0 & x1 != 0 & x01 != 0) {
        f <- "(0 + treatment + treatment:time | cluster)"
    } else if(x0 != 0 & x1 != 0 & x01 == 0) {
        f <- "(0 + treatment + treatment:time || cluster)"
    }
    f
}



# new power
## lme4:::grad.ctr3
gradient <- function (fun, x, delta = 1e-04, ...)
{
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
    if(x0sq == 0 | x1sq == 0) {
        x <- c(x0sq, x1sq)
        x <- x[x > 0]
        x <- sqrt(x)
    } else  {
        x <- matrix(c(x0sq, x01, x01, x1sq), ncol = 2)
        if(x01 == 0) {
            x <- sqrt(x)
            x <- x[c(1,4)]
        } else {
            x <- chol(x)
            x <- x[c(1,3,4)]
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
    ind <- which(para != 0)
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
get_power_new <- function(object, df = "balanced", alpha = 0.05, d = NULL) {

    if(is.null(d)) d <- simulate_data(object)
    f <- lme4::lFormula(formula = create_lmer_formula(object),
                   data = d)

    pc <- setup_power_calc(d, f, object)
    pars <- pc$pars
    X <- pc$X
    Zt <- pc$Zt
    L0 <- pc$L0
    Lambdat <- pc$Lambdat
    Lind <- pc$Lind

    varb <- varb_func(para = pars, X = X, Zt = Zt, L0 = L0, Lambdat = Lambdat, Lind = Lind)
    Phi <- varb(Lc = diag(4))

    if(df == "satterthwaite") {
        df <- get_satterth_df(object, d = d, pars = pars, Lambdat = Lambdat, X = X, Zt = Zt, L0 = L0, Phi = Phi, varb = varb)
    } else if(df == "balanced") {
        df <- get_balanced_df(object)
    } else if(is.numeric(df)) df <- df

    # power

    slope_diff <- get_slope_diff(object)/object$T_end
    se <- sqrt(Phi[4,4])
    lambda <- slope_diff / se

    power <- pt(qt(1-alpha/2, df = df), df = df, ncp = lambda, lower.tail = FALSE) +
        pt(qt(alpha/2, df = df), df = df, ncp = lambda)


    list(power = power, df = df, se = se)

}
