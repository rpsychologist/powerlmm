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
