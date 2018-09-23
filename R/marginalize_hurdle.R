#' @export
marginalize <- function(object, R, ...) {
    UseMethod("marginalize")
}

#' @export
marginalize.plcp_hurdle <- function(object,
                                    R = 1e4,
                                    vectorize = FALSE,
                                    ...) {
    pars <- object

    # pars
    R_cov <- create_R_cov(pars)
    time <- 0:(pars$n1 - 1)

    d <- expand.grid(time = time,
                     treatment = 0:1,
                     subject = 1)

    betas <- with(pars, c(fixed_intercept,
                          fixed_slope,
                          0,
                          log(RR_cont)/pars$T_end))
    betas_hu <- with(pars, c(fixed_hu_intercept,
                             fixed_hu_slope,
                             0,
                             log(OR_hu)/pars$T_end))
    Xmat <- model.matrix(~time * treatment,
                         data = d)


    Zmat <- model.matrix(~time,
                         data = d)

    if(pars$family == "gamma") {
        sd_log <- NULL
        shape <- pars$shape
    } else if(pars$family == "lognormal") {
        sd_log <- pars$sigma_log
        shape <- NULL
    }

    .func <- ifelse(vectorize,
                    ".marginalize_sim_vec",
                    ".marginalize_sim")
    out <- do.call(.func, list(d = d,
                               betas = betas,
                               Xmat = Xmat,
                               Zmat = Zmat,
                               betas_hu = betas_hu,
                               R_cov = R_cov,
                               sd_log = sd_log,
                               shape = shape,
                               marginal = pars$marginal,
                               family = pars$family,
                               R = R,
                               ...))


    out
}

.calc_eta_hurdle <- function(mu, p, marginal, family, sd_log) {
    if(marginal) {
        if(family == "gamma") {
            # Y
            mu_overall <- mu
            # Y > 0
            marg_cont <- mean(exp(mu - log(1 - p)))
            # Median(Y | Y > 0)
            median_cont <- median(exp(mu - log(1 - p)))
        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu
            # Y > 0
            marg_cont <- mean(exp(mu - log(1 - p) - sd_log^2/2))
            # Median(Y | Y > 0)
            median_cont <- median(exp(mu - log(1 - p) - sd_log^2/2))
        }

    }  else {
        if(family == "gamma") {
            # Y
            mu_overall <- mu + log(1 - p)
            # Y > 0
            marg_cont <- mean(exp(mu))
            # Median(Y | Y > 0)
            median_cont <- median(exp(mu))
        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu + log(1 - p) + sd_log^2/2
            # Y > 0
            marg_cont <- mean(exp(mu + sd_log^2/2))
            # Median(Y | Y > 0)
            median_cont <- median(exp(mu + sd_log^2/2))
        }

    }

    list(mu_overall = mu_overall,
         marg_cont = marg_cont,
         median_cont = median_cont)
}

.marginalize_sim <- function(d,
                             betas,
                             betas_hu,
                             R_cov,
                             sd_log,
                             shape,
                             family,
                             marginal = FALSE,
                             R,
                             full = FALSE, ...) {


    #d <- d[sample(1:nrow(d), ceiling(0.7 * nrow(d))), ]

    sd0 <- MASS::mvrnorm(R, c(0,0,0,0), R_cov)
    X <- model.matrix(~time * treatment,
                      data = d)

    Xmat <- X %*% betas
    Xmat_hu <- X %*% betas_hu
    Z <- model.matrix(~time,
                      data = d)

    XtX <- crossprod(X)

    eta_sum <- function(x) {
        cbind("mean" = mean(x),
              "sd" = sd(x),
              "Q2.5" = quantile(x, probs = 0.025),
              "Q25" = quantile(x, probs = 0.25),
              "Q50" = median(x),
              "Q75" = quantile(x, probs = 0.75),
              "Q97.5" = quantile(x, probs = 0.975)
        )

    }
    trans_eta <- function(x, var) {
        out <- do.call(rbind, tmp[, var])
        out <- as.data.frame(out)
        out$time <- d$time
        out$treatment <- d$treatment

        out
    }

    calc_eta <- function(i, full) {

        mu <- Xmat[i, ] + c(Z[i, ] %*% t(sd0[, c(1,2)]))
        hu <- Xmat_hu[i, ] + c(Z[i, ] %*% t(sd0[, c(3,4)]))
        p <- plogis(hu)

        eta <- .calc_eta_hurdle(mu = mu,
                                p = p,
                                marginal = marginal,
                                family = family,
                                sd_log = sd_log)
        mu_overall <- eta$mu_overall
        marg_cont <- eta$marg_cont
        median_cont <- eta$median_cont

        exp_mu_overall <- exp( mu_overall)

            marg_overall <- mean(exp_mu_overall)
            marg_overall_log <- log(marg_overall)
            median_overall <- median(exp_mu_overall)

            marg_cont_log <- log(marg_cont)
            marg_binom <- mean(p)
            marg_binom_logit <- qlogis(marg_binom)


            out <- list(marg_binom = marg_binom,
                        marg_binom_logit = marg_binom_logit,
                        marg_cont= marg_cont,
                        marg_cont_log = marg_cont_log,
                        median_cont = median_cont,
                        marg_overall = marg_overall,
                        marg_overall_log = marg_overall_log,
                        median_overall = median_overall)
        if(full) {
            out <- c(out,
                     list(marg_overall = eta_sum(exp(mu_overall)),
                        p = eta_sum(p),
                        marg_cont = eta_sum(marg_cont))
                     )

        }

        out
    }
    tmp <- lapply(1:nrow(X), calc_eta, full = full)
    tmp <- as.data.frame(do.call(rbind, tmp))

    tmp$time <- d$time
    tmp$treatment <- d$treatment


    # Hedeker et al 2018
    # solve(t(X) %*% X) %*% t(X) %*% tmp$marg_overall

    coefs_median_response <- solve(XtX, crossprod(X, unlist(tmp$median_overall)))
    coefs_median <- solve(XtX, crossprod(X, log(unlist(tmp$median_overall))))
    coef_overall_response <- solve(XtX, crossprod(X, unlist(tmp$marg_overall)))
    coef_overall <- solve(XtX, crossprod(X, unlist(tmp$marg_overall_log)))

    coefs <- cbind(coef_overall)
    colnames(coefs) <- c("overall")

    coefs_median <- cbind(coefs_median)
    colnames(coefs_median) <- c("median_overall")

    coefs <- mapply(function(x, name) {
        x <- data.frame(t(x), check.names = FALSE)
        colnames(x) <- paste(name, names(x), sep = "_")

        x
    },
    list(coef_overall, coefs_median),
    name = c("overall", "median_overall"), SIMPLIFY = FALSE)
    #
    coefs <- do.call(cbind, coefs)
    post <- tmp[tmp$time == max(tmp$time), ]
    marg_post_tx <- mean(post[post$treatment == 1, "marg_overall"])
    marg_post_cc <- mean(post[post$treatment == 0, "marg_overall"])
    median_post_tx <- mean(post[post$treatment == 1, "median_overall"])
    median_post_cc <- mean(post[post$treatment == 0, "median_overall"])
    marg_RR <- marg_post_tx/marg_post_cc
    median_RR <- median_post_tx/median_post_cc

    out <- cbind(coefs,
                 marg_post_tx,
                 marg_post_cc,
                 marg_post_diff = marg_post_tx - marg_post_cc,
                 marg_RR,
                 median_post_tx,
                 median_post_cc,
                 median_post_diff = median_post_tx - median_post_cc,
                 median_RR
    )
    if(full) {

        out <- list(out,
                    #time = tmp,
                    marg_overall = trans_eta(tmp, "marg_overall"),
                    marg_cont = trans_eta(tmp, "marg_cont"),
                    hu_prob = trans_eta(tmp, "p")
                    #coef_overall_response = coef_overall_response,
                    #coefs_median_response = coefs_median_response,
                    #coefs_median = coefs_median,
                    #coef_overall = coef_overall
                    )
    }
    out
}


.marginalize_sim_vec <- function(d,
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
                         c(0,0,0,0),
                         R_cov)

    Xeta <- Xmat %*% betas
    Xeta_hu <- Xmat %*% betas_hu
    mu <- c(Xeta) + tcrossprod(Zmat, sd0[, c(1,2)])
    hu <- c(Xeta_hu) + tcrossprod(Zmat, sd0[, c(3,4)])
    p <- plogis(hu)
    if(marginal) {
        if(family == "gamma") {
            # Y
            mu_overall <- mu
        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu + sd_log^2/2
        }

    }  else {
        if(family == "gamma") {
            # Y
            mu_overall <- mu + log(1 - p)

        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu + log(1 - p) + sd_log^2/2
        }

    }

    exp_mu_overall <- exp( mu_overall)
    d$marg_overall <- matrixStats::rowMeans2(exp_mu_overall)
    d$median_overall <- matrixStats::rowMedians(exp_mu_overall)

    d$marg_p_overall <- matrixStats::rowMeans2(p)
    d$median_p_overall <- matrixStats::rowMedians(p)

    post <- d[d$time == max(d$time),]
    marg_post_tx <- post[post$treatment == 1, "marg_overall"]
    marg_post_cc <- post[post$treatment == 0, "marg_overall"]
    median_post_tx <- post[post$treatment == 1, "median_overall"]
    median_post_cc <-post[post$treatment == 0, "median_overall"]
    marg_RR <- marg_post_tx/marg_post_cc
    median_RR <- median_post_tx/median_post_cc

    # hu
    marg_p_post_tx <- post[post$treatment == 1, "marg_p_overall"]
    marg_p_post_cc <- post[post$treatment == 0, "marg_p_overall"]
    median_p_post_tx <- post[post$treatment == 1, "median_p_overall"]
    median_p_post_cc <-post[post$treatment == 0, "median_p_overall"]
    marg_p_RR <- marg_p_post_tx/marg_p_post_cc
    median_p_RR <- median_p_post_tx/median_p_post_cc
    marg_p_OR <- (marg_p_post_tx/(1-marg_p_post_tx))/(marg_p_post_cc/(1-marg_p_post_cc))
    median_p_OR <- (median_p_post_tx/(1-median_p_post_tx))/(median_p_post_cc/(1-median_p_post_cc))

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

marginalize.brmsfit <- function(object,
                                R = 1e3,
                                .func = .marginalize_sim_vec,
                                ...) {
    fit <- object
    stopifnot(fit$family$family %in% c("hurdle_lognormal", "hurdle_gamma", "custom"))
    family_name <- fit$family$name
    if(!is.null(family_name) && family_name %in% c("hurdle_lognormal_mtp",
                                                   "hurdle_gamma_mtp")) {
        marginal <- TRUE
        family <- ifelse(family_name == "hurdle_lognormal_mtp",
                         "lognormal",
                         "gamma")

    } else {
        family <- ifelse(fit$family$family == "hurdle_lognormal",
                         "lognormal",
                         "gamma")
        marginal <- FALSE
    }

    ss <- posterior_samples(fit, pars = c("^sd",
                                          "^cor"))

    s_d <- brms::standata(fit)
    # pars
    # Continuous
    X <- s_d$X
    beta <- posterior_samples(fit,
                              pars = c("^b_Intercept",
                                       "^b_time",
                                       "^b_treatment",
                                       "^b_time:treatment"))
    names(beta) <- gsub("b_", "", names(beta))
    beta <- beta[colnames(X)]
    stopifnot(length(beta) == ncol(X))

    # Hurdle
    X_hu <- s_d$X_hu
    beta_hu <- posterior_samples(fit,
                                 pars = c("^b_hu_Intercept",
                                          "^b_hu_time",
                                          "^b_hu_treatment",
                                          "^b_hu_time:treatment"))
    names(beta_hu) <- gsub("b_hu_", "", names(beta_hu))
    beta_hu <- beta_hu[colnames(X_hu)]
    stopifnot(length(beta_hu) == ncol(X_hu))

    time <- sort(unique(X[, "time"]))

    d <- expand.grid(time = time,
                     treatment = 0:1,
                     subject = 1)

    if(family == "lognormal") {
        sd_log <- posterior_samples(fit, pars = "^sigma")
    } else if(family == "gamma") {
        shape <- posterior_samples(fit, pars = "^shape")
    }


    marginalize_posterior <- function(j,
                                      R = R,
                                      Xmat,
                                      Zmat,
                                      marginal = marginal,
                                      sd_log = NULL,
                                      shape = NULL,
                                      ...) {

        sd_log <- sd_log[j, ]
        shape <- shape[j, ]
        pars <- list(sd_hu_intercept= ss$sd_subject__hu_Intercept[j],
                     sd_hu_slope = ss$sd_subject__hu_time[j],
                     sd_intercept = ss$sd_subject__Intercept[j],
                     sd_slope = ss$sd_subject__time[j],
                     cor_intercept_slope = ss$cor_subject__Intercept__time[j],
                     cor_intercept_hu_intercept = ss$cor_subject__Intercept__hu_Intercept[j],
                     cor_intercept_hu_slope = ss$cor_subject__Intercept__hu_time[j],
                     cor_slope_hu_intercept = ss$cor_subject__time__hu_Intercept[j],
                     cor_slope_hu_slope = ss$cor_subject__time__hu_time[j],
                     cor_hu_intercept_hu_slope = ss$cor_subject__hu_Intercept__hu_time[j])
        R_cov <- create_R_cov(pars)

        out <- .func(d = d,
                     betas = t(beta[j, ]),
                     betas_hu = t(beta_hu[j, ]),
                     Xmat = Xmat,
                     Zmat = Zmat,
                     R_cov = R_cov,
                     sd_log = sd_log,
                     shape = shape,
                     family = family,
                     marginal = marginal,
                     R = R)

        out

    }

    # we only care about posttest ES
    d <- d[d$time == max(d$time), ]
    Xmat <- model.matrix(~time * treatment,
                         data = d)
    Zmat <- model.matrix(~time,
                         data = d)

    marg <- lapply(1:nrow(ss),
                   FUN = marginalize_posterior,
                   Xmat = Xmat,
                   Zmat = Zmat,
                   R = R,
                   marginal = marginal,
                   family = family)
    marg <- do.call(rbind, marg)
    marg
}
