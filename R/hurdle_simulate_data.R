# TODO
#' add effect sizes to paras
#' add function to plot by tx group




# random effects
create_R_cov <- function(pars) {
    UseMethod("create_R_cov")
}

create_R_cov.plcp_hurdle <- function(pars) {

    # assigns pars to local envir
    for(i in seq_along(pars)) {
        val <- pars[[i]]
        val <- ifelse(is.null(val), 0, val)

        assign(names(pars[i]), val)
    }

    u01 <- sd_intercept * sd_slope * cor_intercept_slope
    u02 <- sd_intercept * sd_hu_intercept * cor_intercept_hu_intercept
    u03 <- sd_intercept * sd_hu_slope * cor_intercept_hu_slope
    u12 <- sd_slope * sd_hu_intercept * cor_slope_hu_intercept
    u13 <- sd_slope * sd_hu_slope * cor_slope_hu_slope
    u23 <- sd_hu_intercept * sd_hu_slope * cor_hu_intercept_hu_slope

    nms <- c("sd_intercept",
             "sd_slope",
             "sd_hu_intercept",
             "sd_hu_slope")

    R_cov <- matrix(c(sd_intercept^2, u01,         u02,               u03,
                      u01,            sd_slope^2,  u12,               u13,
                      u02,            u12,         sd_hu_intercept^2, u23,
                      u03,            u13,         u23,               sd_hu_slope^2),
                    ncol = 4,
                    dimnames = list(nms, nms))

    R_cov
}
sim_hurdle <- function(n1, n2, T_end,
                          fixed_intercept, fixed_hu_intercept,
                          fixed_slope, fixed_hu_slope,
                     sd_hu_intercept, sd_hu_slope,
                     sd_intercept, sd_slope,
                     cor_intercept_slope,
                     cor_intercept_hu_intercept,
                     cor_intercept_hu_slope,
                     cor_slope_hu_intercept,
                     cor_slope_hu_slope,
                     cor_hu_intercept_hu_slope,
                     sigma_log = NULL,
                     shape = NULL,
                     RR_cont = 1,
                     OR_hu = 1,
                     marginal = FALSE,
                     family = "lognormal",
                     ...) {

    time <- seq(0, T_end, length.out = n1)

    stopifnot(family %in% c("lognormal", "gamma"))
    # Random effects vcov matrix
    pars <- list(sd_hu_intercept= sd_hu_intercept,
                 sd_hu_slope = sd_hu_slope,
                 sd_intercept = sd_intercept,
                 sd_slope = sd_slope,
                 cor_intercept_slope = cor_intercept_slope,
                 cor_intercept_hu_intercept = cor_intercept_hu_intercept,
                 cor_intercept_hu_slope = cor_intercept_hu_slope,
                 cor_slope_hu_intercept = cor_slope_hu_intercept,
                 cor_slope_hu_slope = cor_slope_hu_slope,
                 cor_hu_intercept_hu_slope = cor_hu_intercept_hu_slope)

    R_cov <- create_R_cov.plcp_hurdle(pars)
    R <- MASS::mvrnorm(n2,
                       mu = c(0,0,0,0),
                       Sigma = R_cov)

    if(n2 == 1) {
        R <- as.matrix(t(R))
    }

    ## Logistic part
    id <- rep(1:n2, each = length(time))

    b0_hu <- fixed_hu_intercept + R[,3][id]
    b1_hu <- fixed_hu_slope + R[,4][id]
    logit <- b0_hu + b1_hu * time
    yh <- rbinom(n2 * length(time), 1, prob = plogis(logit))

    ## Continuous part
    nh <- which(yh == 0) # 0 is gambling

    b0 <- fixed_intercept + R[,1][id]
    b1 <- fixed_slope + R[,2][id]

    mulog <- b0 + b1 * time

    if(family == "lognormal") {
        stopifnot(!is.null(sigma_log))
        if(marginal) mulog <- mulog - log(1-plogis(logit)) - sigma_log^2/2

        tmp <- rlnorm(n2 * length(time),
                      meanlog = mulog,
                      sdlog = sigma_log)
    } else if(family == "gamma") {
        stopifnot(!is.null(shape))
        if(marginal) mulog <- mulog - log(1-plogis(logit))
        tmp <- rgamma(n2 * length(time),
                      shape = shape,
                      rate = shape / exp(mulog))

    }


    y <- rep(0, length(time) * n2)
    y[nh] <- tmp[nh]

    data.frame(subject = id,
               time = time,
               y = y,
               y_c = y,
               subject_intercept = b0,
               subject_slope = b1,
               subject_intercept_hu = b0_hu,
               subject_slope_hu = b1_hu,
               cluster_intercept = b0,
               cluster_slope = b1,
               cluster_intercept_hu = b0_hu,
               cluster_slope_hu = b1_hu,
               abst = ifelse(yh == 1, 0, 1))
}
sim_hurdle_EMA <- function(n1_obs,
                           n1,
                           n2,
                       fixed_intercept, fixed_hu_intercept,
                       fixed_slope, fixed_hu_slope,
                       sd_hu_intercept, sd_hu_slope,
                       sd_intercept, sd_slope,
                       sd_measures = 0,
                       sd_measures_hu = 0,
                       cor_measures = NULL,
                       cor_intercept_slope = 0,
                       cor_intercept_hu_intercept  = 0,
                       cor_intercept_hu_slope  = 0,
                       cor_slope_hu_intercept  = 0,
                       cor_slope_hu_slope  = 0,
                       cor_hu_intercept_hu_slope  = 0,
                       sigma_log = NULL,
                       shape = NULL,
                       RR_cont = 1,
                       OR_hu = 1,
                       marginal = FALSE,
                       family = "lognormal",
                       ...) {

    time <- rep(0:(n1-1), each = n1_obs)
    time_obs <- seq_along(time)
    stopifnot(family %in% c("lognormal", "gamma"))
    # Random effects vcov matrix
    pars <- list(sd_hu_intercept= sd_hu_intercept,
                 sd_hu_slope = sd_hu_slope,
                 sd_intercept = sd_intercept,
                 sd_slope = sd_slope,
                 cor_intercept_slope = cor_intercept_slope,
                 cor_intercept_hu_intercept = cor_intercept_hu_intercept,
                 cor_intercept_hu_slope = cor_intercept_hu_slope,
                 cor_slope_hu_intercept = cor_slope_hu_intercept,
                 cor_slope_hu_slope = cor_slope_hu_slope,
                 cor_hu_intercept_hu_slope = cor_hu_intercept_hu_slope)
    R_cov <- create_R_cov(pars)
    R_lvl3 <- MASS::mvrnorm(n2,
                       mu = c(0,0,0,0),
                       Sigma = R_cov)
    if(n2 == 1) R_lvl3 <- matrix(R_lvl3, ncol = 4)


    R_cov_lvl2 <- matrix(c(sd_measures^2, sd_measures * sd_measures_hu * cor_measures,
                           sd_measures * sd_measures_hu * cor_measures, sd_measures_hu^2),
                         ncol = 2)
    R_lvl2 <- MASS::mvrnorm(n1 * n2,
                            mu = c(0,0),
                            Sigma =R_cov_lvl2)


    ## Logistic part
    id <- rep(1:n2, each = length(time))

    B0_hu <- rep(fixed_hu_intercept, n2)
    B1_hu <- rep(fixed_hu_slope, n2)
    v0_hu <- R_lvl3[,3]
    u0_hu <- rep(R_lvl2[, 2], each = n1_obs)
    v1_hu <- R_lvl3[,4][id]
    logit <- B0_hu[id] + v0_hu[id] + u0_hu + (B1_hu[id] + v1_hu) * time
    yh <- rbinom(n2 * length(time), 1, prob=1/(1 + exp(logit)))

    ## lognormal
    nh <- which(yh == 1) # 0 is gambling

    b0 <- rep(fixed_intercept, n2)
    b1 <- rep(fixed_slope, n2)

    v0 <- R_lvl3[, 1][id]
    v1 <- R_lvl3[, 2][id]
    u0 <- rep(R_lvl2[, 1], each = n1_obs)
    mulog <- b0[id] + (b1[id] + v1) * time + v0 + u0



    if(family == "lognormal") {
        stopifnot(!is.null(sigma_log))
        if(marginal) mulog <- mulog - log(1-plogis(logit)) - sigma_log^2/2

        tmp <- rlnorm(n2 * length(time),
                      meanlog = mulog,
                      sdlog = sigma_log)
    } else if(family == "gamma") {
        stopifnot(!is.null(shape))
        if(marginal) mulog <- mulog - log(1-plogis(logit))
        tmp <- rgamma(n2 * length(time),
                      shape = shape,
                      rate = shape / exp(mulog))

    }


    y <- rep(0, length(time) * n2)
    y[nh] <- tmp[nh]

    data.frame(subject = id,
               time = time,
               time_obs = time_obs,
               mu_log = mulog,
               mu_fixed = b0[id] + (b1[id]) * (time_obs-1)/n1_obs,
               mu_subjects = b0[id] + (b1[id] + v1) * (time_obs-1)/n1_obs + v0,
               hu_logit = logit,
               y = y,
               abst = ifelse(yh == 1, 0, 1))


}
#' Generate Hurdle-lognormal data
#'
#' @param paras paras
#' @param gen_fake func
#'
#' @return df
#' @export
simulate_data.plcp_hurdle <- function(paras,
                                      gen_fake = "sim_hurdle") {
    if (is.data.frame(paras))
        paras <- as.list(paras)
    if(is.null(paras$prepared)) {
        tmp <- prepare_paras(paras)
    } else tmp <- paras
    paras <- tmp$control
    paras_tx <- tmp$treatment

    if(paras$EMA) gen_fake <- "sim_hurdle_EMA"

    d_c <- do.call(gen_fake, paras)
    d_c$treatment <- 0

    paras_tx$fixed_slope <- paras$fixed_slope + log(paras$RR_cont)/paras$T_end
    paras_tx$fixed_hu_slope <- paras$fixed_hu_slope + log(paras$OR_hu)/paras$T_end

    d_tx <-  do.call(gen_fake, paras_tx)
    d_tx$treatment <- 1
    d_tx$subject = d_tx$subject + paras$n2



    # missing data
    if (is.list(paras$dropout) |
        is.function(paras$dropout) |
        length(paras$dropout) > 1) {

        miss_c <- create_dropout_indicator(paras)
        d_c <- add_NA_values_from_indicator(d_c, miss_c)
    }
    if (is.list(paras_tx$dropout) |
        is.function(paras_tx$dropout) |
        length(paras_tx$dropout) > 1) {
        miss_tx <- create_dropout_indicator(paras_tx)
        d_tx <- add_NA_values_from_indicator(d_tx, miss_tx)
    }
    d <- rbind(d_c, d_tx)
    #d$treatment <- factor(d$treatment, labels = c("Control", "Treatment"))

    d
}







