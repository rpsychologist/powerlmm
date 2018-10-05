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
            mu_positive <- exp(mu - log(1 - p))
        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu
            # Y > 0
            marg_cont <- exp(mu - log(1 - p) - sd_log^2/2)
        }

    }  else {
        if(family == "gamma") {
            # Y
            mu_overall <- mu + log(1 - p)
            # Y > 0
            mu_positive <- exp(mu)
        } else if(family == "lognormal") {
            # Y
            mu_overall <- mu + log(1 - p) + sd_log^2/2
            # Y > 0
            mu_positive <- exp(mu + sd_log^2/2)
        }

    }

    list(mu_overall = mu_overall,
         mu_positive = mu_positive)
}

# return summaries of marginal ests distribution
eta_sum <- function(x) {
    cbind("mean" = mean(x),
          "sd" = sd(x),
          "Q0.5" = quantile(x, probs = 0.005),
          "Q2.5" = quantile(x, probs = 0.025),
          "Q10" = quantile(x, probs = 0.10),
          "Q25" = quantile(x, probs = 0.25),
          "Q50" = median(x),
          "Q75" = quantile(x, probs = 0.75),
          "Q90" = quantile(x, probs = 0.90),
          "Q97.5" = quantile(x, probs = 0.975),
          "Q99.5" = quantile(x, probs = 0.995)

    )

}

# helper to display marginal ests
trans_eta <- function(x, var, d) {
    out <- do.call(rbind, x[, var])
    out <- as.data.frame(out)
    out <- cbind(data.frame(var = var,
                            treatment = d$treatment,
                            time = d$time),
                 out)
    rownames(out) <- NULL
    out
}

trans_post_ps <- function(x, hu = FALSE) {
    post_ps <- x[vapply(x, is.data.frame, logical(1))]

    ind <- vapply(post_ps, function(d) all(d$treatment == 0), logical(1))
    post_ps_c <- post_ps[[which(ind)]]
    post_ps_tx <- post_ps[[which(!ind)]]

    if(hu) {
        post_ps <- data.frame(percentile = post_ps_c$percentile,
                              diff = post_ps_tx$value - post_ps_c$value,
                              OR = get_OR(post_ps_tx$value, post_ps_c$value)
                              )
    } else {
        post_ps <- data.frame(percentile = post_ps_c$percentile,
                              diff = post_ps_tx$value - post_ps_c$value,
                              ratio = post_ps_tx$value / post_ps_c$value)
    }


    list(control = post_ps_c,
         treatment = post_ps_tx,
         effect = post_ps)
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
        mu_positive <- eta$mu_positive

        exp_mu_overall <- exp( mu_overall)

        if(i %in% which(d$time == max(d$time))) {
            ps <- 1:99/100
            post <- data.frame("percentile" = ps,
                               "value" = quantile(exp_mu_overall, ps),
                               "treatment" = d[i, "treatment"]
                               )
            post_hu <- data.frame("percentile" = ps,
                               "value" = quantile(p, ps),
                               "treatment" = d[i, "treatment"]
            )

        } else {
            post <- NULL
            post_hu <- NULL
            }

        out <- list(hu_prob = eta_sum(p),
                    marg_y_positive = eta_sum(mu_positive),
                    marg_y_overall = eta_sum(exp_mu_overall),
                    post = post,
                    post_hu = post_hu,
                    exp_mu_overall_vec = exp_mu_overall)


        out
    }
    tmp <- lapply(1:nrow(X), calc_eta, full = full)
    tmp <- as.data.frame(do.call(rbind, tmp))

    post_ps <- trans_post_ps(tmp$post)
    post_hu_ps <- trans_post_ps(tmp$post_hu, hu = TRUE)

    marg_y_overall <- trans_eta(tmp, "marg_y_overall", d = d)
    marg_y_positive <- trans_eta(tmp, "marg_y_positive", d= d)
    hu_prob <- trans_eta(tmp, "hu_prob", d = d)

    # Hedeker et al 2018
    # solve(t(X) %*% X) %*% t(X) %*% tmp$marg_overall
    coef_overall_median_log <- solve(XtX, crossprod(X, log(marg_y_overall[, "Q50"] )))
    coef_overall_marg_log <- solve(XtX, crossprod(X, log(marg_y_overall[, "mean"])))
    coef_hu_prob_marg_logit <- solve(XtX, crossprod(X, qlogis(hu_prob[, "mean"])))
    coef_hu_prob_median_logit <- solve(XtX, crossprod(X, qlogis(hu_prob[, "Q50"])))

    coefs <- mapply(function(x, name) {
            x <- x
            d <- data.frame(var = paste("b", name, rownames(x), sep = "_"),
                            est = c(x),
                            check.names = FALSE)
            d
        },
        list(coef_overall_marg_log,
             coef_overall_median_log,
             coef_hu_prob_marg_logit,
             coef_hu_prob_median_logit),
        name = c("overall_marg",
                 "overall_median",
                 "hu_prob_marg",
                 "hu_prob_median"),
        SIMPLIFY = FALSE)
    #
    names(coefs) <- c("marginal",
                      "median",
                      "marginal",
                      "median")

    coefs <- list("y_overall" = coefs[1:2],
                  "hu_prob" = coefs[3:4])

    ## TODO: also return post ES with sd and percentiles

    # posttest
    post <- marg_y_overall[marg_y_overall$time == max(marg_y_overall$time), c("treatment","mean", "Q50")]
    marg_post_tx <- post[post$treatment == 1, "mean"]
    marg_post_cc <- post[post$treatment == 0, "mean"]
    median_post_tx <- post[post$treatment == 1, "Q50"]
    median_post_cc <- post[post$treatment == 0, "Q50"]
    marg_RR <- marg_post_tx/marg_post_cc
    median_RR <- median_post_tx/median_post_cc

    post <- rbind(marg_post_tx,
                  marg_post_cc,
                  marg_post_diff = marg_post_tx - marg_post_cc,
                  marg_RR,
                  median_post_tx,
                  median_post_cc,
                  median_post_diff = median_post_tx - median_post_cc,
                  median_RR
                  )
    post <- data.frame(var = rownames(post),
                       est = post, row.names = NULL)


    ## post hu, y = 0
    post_hu <- hu_prob[hu_prob$time == max(hu_prob$time), c("treatment","mean", "Q50")]
    marg_hu_post_tx <- post_hu[post_hu$treatment == 1, "mean"]
    marg_hu_post_cc <- post_hu[post_hu$treatment == 0, "mean"]
    median_hu_post_tx <- post_hu[post_hu$treatment == 1, "Q50"]
    median_hu_post_cc <- post_hu[post_hu$treatment == 0, "Q50"]
    marg_OR <- get_OR(marg_hu_post_tx, marg_hu_post_cc)
    median_OR <- get_OR(median_hu_post_tx, median_hu_post_cc)

    post_hu <- rbind(marg_hu_post_tx,
                  marg_hu_post_cc,
                  marg_hu_post_diff = marg_hu_post_tx - marg_hu_post_cc,
                  marg_OR,
                  median_hu_post_tx,
                  median_hu_post_cc,
                  median_hu_post_diff = median_hu_post_tx - median_hu_post_cc,
                  median_OR
    )
    post_hu <- data.frame(var = rownames(post_hu),
                          est = post_hu,
                          row.names = NULL)

    list(coefs = coefs,
         y_overall = marg_y_overall,
         y_positive = marg_y_positive,
         hu_prob = hu_prob,
         post = post,
         post_hu = post_hu,
         post_ps = post_ps,
         post_hu_ps = post_hu_ps,
         mu_overall_vec = tmp$exp_mu_overall_vec
         )


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
        R_cov <- create_R_cov.plcp_hurdle(pars)

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
