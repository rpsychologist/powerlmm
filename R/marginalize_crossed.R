#' @export
marginalize.plcp_crossed <- function(object,
                                    R = 1e4,
                                    vectorize = FALSE,
                                    ...) {
    pars <- object
    ## Vcov level-2
    R_cov2 <- get_lvl2_vcov(pars)
    R_cov2[is.na(R_cov2)] <- 0
    ## Vcov level-3
    R_cov3 <- get_lvl3_vcov(pars)
    R_cov3[is.na(R_cov3)] <- 0
    d <- .create_dummy_d(pars)
    betas <- with(pars, c(fixed_intercept,
                          fixed_slope,
                          0,
                          get_slope_diff(pars)/pars$T_end))
    X <- model.matrix(~time * treatment,
                         data = d)
    Z2 <- model.matrix(~time,
                         data = d)
    Z3 <- X
    if(pars$family == "gamma") {
        sigma_error <- NULL
        shape <- pars$shape
    } else {
        sigma_error <- pars$sigma_error
        shape <- NULL
    }
    .func <- ifelse(vectorize,
                    ".marginalize_nested_sim_vec",
                    ".marginalize_nested_sim")
    out <- do.call(.func, list(d = d,
                               betas = betas,
                               X = X,
                               Z2 = Z2,
                               Z3 = Z3,
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