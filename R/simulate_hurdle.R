# Methods for simulating a hurdle model

get_balanced_df.plcp_hurdle <- function(object) {
    100
}
#' #' @export
#' get_tot_n.plcp_hurdle <- function(paras) {
#'     paras$n2 * 2
#' }
#' @export
get_RE_thetas.plcp_hurdle <- function(paras) {
    data.frame(
        parameter = c(
            "subject_intercept",
            "subject_slope",
            "cluster_intercept",
            "cluster_intercept_tx",
            "cluster_slope",
            "cluster_slope_tx",
            "error",
            "cor_subject",
            "cor_cluster_intercept_slope",
            "cor_cluster_intercept_intercept_tx",
            "cor_cluster_intercept_slope_tx",
            "cor_cluster_intercept_tx_slope_tx",
            "cor_cluster_slope_intercept_tx",
            "cor_cluster_slope_slope_tx"
        ),
        theta = 0
    )

}
#' @export
get_slope_diff.plcp_hurdle <- function(paras, hu = FALSE) {
    if(hu) {
        log(paras$OR_hu)
    } else {
        log(paras$RR_cont)
    }
}

#' #' @export
#' get_n3.plcp_hurdle <- function(paras) {
#'     NA
#' }
#' @export
# get_n2.plcp_hurdle <- function(paras) {
#     treatment <- 0
#     attr(treatment, "per_treatment") <- FALSE
#     control <- 0
#     attr(control, "per_treatment") <- FALSE
#     list(treatment = treatment,
#          control = control)
# }
#' @export
get_FE_thetas.plcp_hurdle <- function(paras, marginalize = FALSE ,... ) {
    out <- list(
        "Intercept" = paras$fixed_intercept,
        "treatment" = 0,
        "time" = paras$fixed_slope,
        "time:treatment" = get_slope_diff(paras)/paras$T_end,
        "hu_Intercept" = paras$fixed_hu_intercept,
        "hu_treatment" = 0,
        "hu_time" = paras$fixed_hu_slope,
        "hu_time:treatment" = get_slope_diff(paras, hu = TRUE)/paras$T_end
    )

    if(marginalize) {
        out <- c(out,
                 as.list( marginalize(paras, R = 1e5))
        )
    }

    out
}

#' @export
get_RE_thetas.plcp_hurdle <- function(paras) {
    list(
        "sd_subject__Intercept" = paras$sd_intercept,
        "sd_subject__hu_Intercept" = paras$sd_hu_intercept,
        "sd_subject__time" = paras$sd_slope,
        "sd_subject__hu_time" = paras$sd_hu_slope,
        "cor_subject__Intercept__time" = paras$cor_intercept_slope,
        "cor_subject__Intercept__hu_Intercept" = paras$cor_intercept_hu_intercept,
        "cor_subject__Intercept__hu_time" = paras$cor_intercept_hu_slope,
        "cor_subject__time__hu_Intercept" = paras$cor_slope_hu_intercept,
        "cor_subject__time__hu_time" = paras$cor_slope_hu_slope,
        "cor_subject__hu_Intercept__hu_time" = paras$cor_hu_intercept_hu_slope
    )
}
#' @export
sim_formula.brmsfit <- function(formula,
                                data_transform = NULL,
                                test = "time:treatment",
                                marginalize = FALSE,
                                chains = 1,
                                iter = 2000,
                                refresh = 0,
                                silent = TRUE,
                                cores = 1,
                                ...) {
    x <- list("formula" = formula,
              "data_transform" = data_transform,
              "data_transform_lab" = substitute(data_transform),
              "test" = test,
              "chains" = chains,
              "iter" = iter,
              "refresh" = refresh,
              "silent" = silent,
              "cores" = 1,
              "marginalize" = marginalize,
              ...)

    class(x) <- append(class(x), c("plcp_brmsformula", "plcp_sim_formula"))
    x
}

#' @export
fit_model.plcp_brmsformula <- function(formula,
                                       data,
                                       ...) {
    args <- list(object = formula$formula,
                newdata = data,
                chains = formula$chains,
                iter = formula$iter,
                refresh = formula$refresh,
                cores = formula$cores)
    if(formula$silent) {
        capture.output(out <- do.call(update, args))
    } else {
        out <- do.call(update, args)
    }

    out
}
#' @export
get_fixef_coef.brmsfit <- function(fit, formula, ...) {

    out <- brms::fixef(fit, robust = TRUE)[, c("Estimate")]

    if(formula$marginalize) {
        marg <- marginalize(fit, R = 1e3)
        marg <- apply(marg, 2, median)
        out <- c(out, marg)
    }
    out
}
#' @export
get_fixef.brmsfit <- function(fit, test, df_bw, satterthwaite, formula) {

    # need to know in which order time and treatment was entered
    # then adjust 'fit$test' to match
    #FE_coefs <- get_fixef_coef(fit, ...)

    FE_coefs <- brms::fixef(fit)
    colnames(FE_coefs) <- c("estimate",
                            "se",
                            "CI_lwr",
                            "CI_upr")
    FE_coefs <- cbind(data.frame(parameter = rownames(FE_coefs),
                                 stringsAsFactors = FALSE),
                      FE_coefs)

    if(formula$marginalize) {
        marg <- marginalize(fit, R = 1e3)
        marg <- data.frame("parameter" = colnames(marg),
                           "estimate" = apply(marg, 2, median),
                           "se" = apply(marg, 2, sd),
                           "CI_lwr" = apply(marg, 2, quantile, 0.025),
                           "CI_upr" = apply(marg, 2, quantile, 0.975)
        )
        FE_coefs <- rbind(FE_coefs, marg)

    }

    rownames(FE_coefs) <- NULL

    class(FE_coefs) <- append("plcp_brmsformula", class(FE_coefs))
    FE_coefs
}


#' @export
add_p_value.brmsfit <- function(fit, ...) {
    list("df" = NA,
         "p" = NA)
}

#' @export
extract_random_effects.brmsfit <- function(fit) {
    ss <- brms::posterior_samples(fit, pars =  c("^sd",
                                           "^cor",
                                           "^sigma",
                                           "^shape"))
    data.frame(parameter = colnames(ss),
               var1 = NA,
               var2 = NA,
               vcov = apply(ss, 2, mean),
               CI_lwr = apply(ss, 2, quantile, probs = 0.025),
               CI_upr = apply(ss, 2, quantile, probs = 0.975),
               stringsAsFactors = FALSE)

}
#' @export
get_convergence.brmsfit <- function(fit) {
    np <- brms::nuts_params(fit)

    list("divergent" = sum(subset(np, Parameter == "divergent__")$Value))
}
#' @export
rename_random_effects.plcp_brmsformula <- function(.x, ...) {
    # currently not renamed
    .x
}

#' @export
get_CI.brmsfit <- function(fit, test, FE, ...) {
    # CI included by default; get_fixef.brmsfit()
    FE
}

#' @export
get_LL.brmsfit <- function(fit) {
    list(ll = NA,
         df = NA)
}

