# some code for the custom brms marginal hurdle models

#' @export
hurdle_gamma_mtp <- brms::custom_family(
    "hurdle_gamma_mtp",
    dpars = c("shape","mu", "hu"),
    lb = c(0, NA, 0),
    ub = c(NA, NA, 1),
    links = c("log", "identity", "logit"),
    type = "real"
)


# Adapted from brms hurdle_gamma family
# create by Paul-Christian Buerkner
hurdle_gamma_mtp_scode <- "
/* marginal hurdle gamma log-PDF of a single response
* logit parameterization of the hurdle part
* Args:
*   y: the response value
*   alpha: shape parameter of the gamma distribution
*   mu: linear predictor for the overall mean
*   hu: linear predictor for the hurdle part
* Returns:
*   a scalar to be added to the log posterior
*/
real hurdle_gamma_mtp_lpdf(real y, real alpha, real mu, real hu) {
if (y == 0) {
return bernoulli_logit_lpmf(1 | hu);
} else {
real beta;
beta = alpha * exp(-(mu - log1m_inv_logit(hu)));
return bernoulli_logit_lpmf(0 | hu) +
gamma_lpdf(y | alpha, beta);
}
}
"



#' @export
hurdle_gamma_mtp_stanvars <- brms::stanvar(scode = hurdle_gamma_mtp_scode,
                                           block = "functions")

#' @export
predict_hurdle_gamma_mtp <- function (i, draws, ...) {

    # Theta proportions
    theta <- brms:::get_dpar(draws, "hu", i = i)
    shape <- brms:::get_dpar(draws, "shape", i = i)

    # mu on log scale?
    # solve for g(Y | y > 0)
    mu <- brms:::get_dpar(draws, "mu", i = i) - log(1-theta)

    ndraws <- draws$nsamples
    hu <- runif(ndraws, 0, 1)
    ifelse(hu < theta, 0, rgamma(ndraws, shape = shape, rate = shape / exp(mu) ))
}

#' @export
fitted_hurdle_gamma_mtp <- function (draws, ...) {

    with(draws$dpars, exp(mu))
}
