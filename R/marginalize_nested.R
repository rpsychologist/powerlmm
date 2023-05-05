# Marginalize nested GLMMs

# TODO
# * use correct inv_link
# * update sim_vec version

# @param p1 prob tx
# @param p0 prob control
get_OR <- function(p1, p0) {
    (p1 / (1 - p1)) / (p0 / (1 - p0))
}

#' @export
#' @importFrom stats quantile plogis qlogis model.matrix
#' rnorm rlogis rpois rgamma rlnorm rbinom
marginalize.plcp_nested <- function(object,
                                    R = 1e4,
                                    vectorize = FALSE,
                                    ...) {
    pars <- object

    # pars
    ## Vcov level-2
    Sigma_subject <- with(pars,
        c(
            sigma_subject_intercept^2,
            sigma_subject_intercept * sigma_subject_slope * cor_subject,
            sigma_subject_intercept * sigma_subject_slope * cor_subject,
            sigma_subject_slope^2
        )
    )
    R_cov2 <- matrix(Sigma_subject, 2, 2)
    R_cov2[is.na(R_cov2)] <- 0
    ## Vcov level-3

    Sigma_cluster <- with(pars,
        c(
            sigma_cluster_intercept^2,
            sigma_cluster_intercept * sigma_cluster_slope * cor_cluster,
            sigma_cluster_intercept * sigma_cluster_slope * cor_cluster,
            sigma_cluster_slope^2
        )
    )
    R_cov3 <- matrix(Sigma_cluster, 2, 2)
    R_cov3[is.na(R_cov3)] <- 0
    d <- .create_dummy_d(pars)
    betas <- with(pars, c(fixed_intercept,
        fixed_slope,
        0,
        get_slope_diff(pars) / pars$T_end))
    X <- model.matrix(~ time * treatment,
        data = d)
    Z <- model.matrix(~time,
        data = d)
    sigma_error <- pars$sigma_error
    shape <- NULL
    .func <- ifelse(vectorize,
        ".marginalize_nested_sim_vec",
        ".marginalize_nested_sim")
    out <- do.call(.func, list(d = d,
        betas = betas,
        X = X,
        Z2 = Z,
        Z3 = Z,
        R_cov2 = R_cov2,
        R_cov3 = R_cov3,
        sigma = sigma_error,
        shape = shape,
        family = pars$family,
        partially_nested = pars$partially_nested,
        R = R,
        ...))

    out$paras <- object
    class(out) <- c("plcp_marginal_nested")
    out
}

.create_dummy_d <- function(pars) {
    time <- get_time_vector(pars)

    d <- expand.grid(time = time,
        treatment = 0:1,
        subject = 1)

    d
}

.get_inv_link <- function(family = "gaussian", sigma = NULL) {
    # log-transformation, but call link for concistency
    ln_inv <- function(sigma) {
        sigma <- sigma
        function(eta) {
            exp(eta + sigma^2 / 2)
        }
    }

    inv_link <- switch(as.character(family),
        "gaussian" = function(eta) eta,
        "binomial" = plogis,
        "poisson" = exp,
        "gamma" = exp,
        "lognormal" = ln_inv(sigma)
    )

    inv_link
}
.get_link <- function(family, sigma = NULL) {
    ln_link <- function(sigma) {
        sigma <- sigma
        function(y) {
            log(y) - sigma^2 / 2
        }
    }


    link <- switch(as.character(family),
        "gaussian" = function(eta) eta,
        "binomial" = qlogis,
        "poisson" = log,
        "gamma" = log,
        "lognormal" = ln_link(sigma)
    )

    link
}


# marginalize ests over random effects
.marginalize_nested_sim <- function(d,
                                    betas,
                                    X,
                                    Z2,
                                    Z3,
                                    R_cov2,
                                    R_cov3,
                                    sigma,
                                    shape,
                                    partially_nested,
                                    R,
                                    link_scale = FALSE,
                                    full = FALSE,
                                    ...) {
    sd2 <- MASS::mvrnorm(R, rep(0, ncol(R_cov2)), R_cov2)
    sd3 <- MASS::mvrnorm(R, rep(0, ncol(R_cov3)), R_cov3)
    Xmat <- X %*% betas
    XtX <- crossprod(X)
    family <- "gaussian"
    inv_link <- .get_inv_link(family, sigma)
    link <- .get_link(family, sigma)
    if (link_scale) {
        inv_link <- function(eta) eta
        link <- function(eta) eta
    }
    calc_eta <- function(i) {
        if (partially_nested) {
            tx <- d[i, "treatment"]
            sd3 <- sd3 * tx # cc = 0
            sd0 <- sd2 + sd3
        }
        # level 2 (includes lvl 3)
        v <- c(Z3[i, ] %*% t(sd3))
        mu2 <- Xmat[i, ] + c(Z2[i, ] %*% t(sd2)) + v
        # level 3 only
        mu3 <- Xmat[i, ] + v
        exp_mu2 <- inv_link(mu2)
        exp_mu3 <- inv_link(mu3)
        if (i %in% which(d$time == max(d$time))) {
            ps <- (1:99) / 100
            post <- data.frame("percentile" = ps,
                "value" = quantile(exp_mu2, ps),
                "treatment" = d[i, "treatment"]
            )
        } else {
            post <- NULL
        }
        out <- list(marg_y2 = eta_sum(exp_mu2),
            marg_y3 = eta_sum(exp_mu3),
            post = post,
            exp_mu2_vec = exp_mu2,
            exp_mu3_vec = exp_mu3)

        out
    }
    tmp <- lapply(1:nrow(X), calc_eta)
    tmp <- as.data.frame(do.call(rbind, tmp))

    post_ps <- trans_post_ps(tmp$post, hu = family == "binomial")
    # post_hu_ps <- trans_post_ps(tmp$post_hu, hu = TRUE)

    marg_y <- trans_eta(tmp, "marg_y2", d = d)
    marg_y3 <- trans_eta(tmp, "marg_y3", d = d)

    # Hedeker et al 2018
    # solve(t(X) %*% X) %*% t(X) %*% tmp$marg_overall
    coef_median_log <- solve(XtX, crossprod(X, link(marg_y[, "Q50"])))
    coef_marg_log <- solve(XtX, crossprod(X, link(marg_y[, "mean"])))

    coefs <- mapply(function(x, name) {
        x <- x
        d <- data.frame(var = paste("b", name, rownames(x), sep = "_"),
            est = c(x),
            check.names = FALSE)
        d
    },
    list(coef_median_log,
        coef_marg_log),
    name = c("y_median",
        "y_marg"),
    SIMPLIFY = FALSE)
    #
    names(coefs) <- c("median",
        "marginal")
    coefs <- list("y" = coefs)

    # posttest
    post <- marg_y[marg_y$time == max(marg_y$time), c("treatment", "mean", "Q50")]
    marg_post_tx <- post[post$treatment == 1, "mean"]
    marg_post_cc <- post[post$treatment == 0, "mean"]
    median_post_tx <- post[post$treatment == 1, "Q50"]
    median_post_cc <- post[post$treatment == 0, "Q50"]
    marg_RR <- marg_post_tx / marg_post_cc
    median_RR <- median_post_tx / median_post_cc
    marg_OR <- get_OR(marg_post_tx, marg_post_cc)
    median_OR <- get_OR(median_post_tx, median_post_cc)


    post <- rbind(marg_post_tx,
        marg_post_cc,
        marg_post_diff = marg_post_tx - marg_post_cc,
        marg_RR,
        marg_OR,
        median_post_tx,
        median_post_cc,
        median_post_diff = median_post_tx - median_post_cc,
        median_RR,
        median_OR
    )
    post <- data.frame(var = rownames(post),
        est = post, row.names = NULL)
    list(coefs = coefs,
        y2 = marg_y,
        y3 = marg_y3,
        post = post,
        post_ps = post_ps,
        mu2_vec = tmp$exp_mu2_vec,
        mu3_vec = tmp$exp_mu3_vec
    )
}


.marginalize_nested_sim_vec <- function(d,
                                        betas,
                                        betas_hu,
                                        R_cov,
                                        sd_log,
                                        shape,
                                        family,
                                        Xmat,
                                        Zmat,
                                        marginal = FALSE,
                                        R,
                                        full = FALSE, ...) {


    sd0 <- MASS::mvrnorm(R,
        c(0, 0, 0, 0),
        R_cov)

    Xeta <- Xmat %*% betas
    Xeta_hu <- Xmat %*% betas_hu
    mu <- c(Xeta) + tcrossprod(Zmat, sd0[, c(1, 2)])
    hu <- c(Xeta_hu) + tcrossprod(Zmat, sd0[, c(3, 4)])
    p <- plogis(hu)
    if (marginal) {
        if (family == "gamma") {
            # Y
            mu_overall <- mu
        } else if (family == "lognormal") {
            # Y
            mu_overall <- mu + sd_log^2 / 2
        }
    } else {
        if (family == "gamma") {
            # Y
            mu_overall <- mu + log(1 - p)
        } else if (family == "lognormal") {
            # Y
            mu_overall <- mu + log(1 - p) + sd_log^2 / 2
        }
    }
    exp_mu_overall <- exp(mu_overall)
    d$marg_overall <- matrixStats::rowMeans2(exp_mu_overall)
    d$median_overall <- matrixStats::rowMedians(exp_mu_overall)

    d$marg_p_overall <- matrixStats::rowMeans2(p)
    d$median_p_overall <- matrixStats::rowMedians(p)

    post <- d[d$time == max(d$time), ]
    marg_post_tx <- post[post$treatment == 1, "marg_overall"]
    marg_post_cc <- post[post$treatment == 0, "marg_overall"]
    median_post_tx <- post[post$treatment == 1, "median_overall"]
    median_post_cc <- post[post$treatment == 0, "median_overall"]
    marg_RR <- marg_post_tx / marg_post_cc
    median_RR <- median_post_tx / median_post_cc

    # hu
    marg_p_post_tx <- post[post$treatment == 1, "marg_p_overall"]
    marg_p_post_cc <- post[post$treatment == 0, "marg_p_overall"]
    median_p_post_tx <- post[post$treatment == 1, "median_p_overall"]
    median_p_post_cc <- post[post$treatment == 0, "median_p_overall"]
    marg_p_RR <- marg_p_post_tx / marg_p_post_cc
    median_p_RR <- median_p_post_tx / median_p_post_cc
    marg_p_OR <- get_OR(marg_p_post_tx, marg_p_post_cc)
    median_p_OR <- get_OR(median_p_post_tx, median_p_post_cc)

    out <- cbind(marg_post_tx,
        marg_post_cc,
        marg_post_diff = marg_post_tx - marg_post_cc,
        marg_RR,
        median_post_tx,
        median_post_cc,
        median_post_diff = median_post_tx - median_post_cc,
        median_RR,
        marg_p_post_tx,
        marg_p_post_cc,
        marg_p_post_diff = marg_p_post_tx - marg_p_post_cc,
        marg_p_RR,
        marg_p_OR,
        median_p_post_tx,
        median_p_post_cc,
        median_p_post_diff = median_p_post_tx - median_p_post_cc,
        median_p_RR,
        median_p_OR
    )
    out
}


## Sample level 1
.sample_level1_nested <- function(pars,
                                  fixed_subject_percentiles = c(0.5, 0.5),
                                  fixed_cluster_percentiles = c(0.5, 0.5),
                                  R = 1e4,
                                  link_scale = FALSE,
                                  ...) {

    # arsg:
    #   fixed_subject_percentiles: level 2 percentiles
    #   fixed_cluster_percentiles: level 2 percentiles

    d <- .create_dummy_d(pars)
    family <- pars$family

    betas <- with(pars, c(fixed_intercept,
        fixed_slope,
        0,
        get_slope_diff(pars) / pars$T_end))

    X <- model.matrix(~ time * treatment,
        data = d)

    Xmat <- X %*% betas
    Z <- model.matrix(~time,
        data = d)

    XtX <- crossprod(X)

    inv_link <- .get_inv_link(family, pars$sigma_error)
    link <- .get_link(family, pars$sigma_error)

    # if(link_scale) {
    #    inv_link <- function(eta) eta
    #    link <- function(eta) eta
    # }

    # RE
    sd3 <- qnorm(fixed_cluster_percentiles, 0, with(pars, c(sigma_cluster_intercept,
        sigma_cluster_slope)))
    sd2 <- qnorm(fixed_subject_percentiles,
        0,
        with(pars, c(sigma_subject_intercept,
            sigma_subject_slope)))

    sd3[is.na(sd3)] <- 0
    sd2[is.na(sd2)] <- 0
    sd0 <- sd2 + sd3

    calc_eta <- function(i, full) {
        if (pars$partially_nested) {
            tx <- d[i, "treatment"]
            sd3 <- sd3 * tx # cc = 0
            sd0 <- sd2 + sd3
        }

        mu <- Xmat[i, ] + Z[i, ] %*% sd0

        inv_mu <- inv_link(mu)
        if (family == "gaussian") {
            y <- rnorm(R, inv_mu, pars$sigma_error)
        } else if (family == "binomial") {
            y <- rlogis(R, location = mu)
            if (!link_scale) y <- as.numeric(y > 0)
        } else if (family == "poisson") {
            y <- rpois(R, lambda = inv_mu)
            if (link_scale) y <- log(y)
        }
        else if (family == "gamma") {
            shape <- pars$shape
            y <- rgamma(R,
                shape = shape,
                rate = shape / inv_mu)
            if (link_scale) y <- log(y)
        } else if (family == "lognormal") {
            y <- rlnorm(R,
                meanlog = mu,
                sdlog = pars$sigma_error)
            if (link_scale) y <- log(y)
        }

        out <- list(marg_y1 = eta_sum(y),
            exp_mu1_vec = y)
    }
    tmp <- lapply(1:nrow(X), calc_eta)
    tmp <- as.data.frame(do.call(rbind, tmp))

    marg_y <- trans_eta(tmp, "marg_y1", d = d)

    list(coefs = NULL,
        y = marg_y,
        post = NULL,
        post_ps = NULL,
        mu1_vec = tmp$exp_mu1_vec
    )

}