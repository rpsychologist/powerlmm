#' @method study_parameters plcp_design_hurdle_nested
#' @export
study_parameters.plcp_design_hurdle_nested <- function(design = NULL,
                                         n1_obs = NULL, # obs per time point (EMA)
                                         n1 = NULL, # time points
                                         n2,
                                         T_end = NULL,
                                         fixed_intercept,
                                         fixed_intercept_tx = 0,
                                         fixed_hu_intercept,
                                         fixed_hu_intercept_tx = 0,
                                         fixed_slope,
                                         fixed_hu_slope,
                                         sd_hu_intercept = 0,
                                         sd_hu_slope = 0,
                                         sd_intercept = 0,
                                         sd_slope = 0,
                                         cor_intercept_slope = 0,
                                         cor_intercept_hu_intercept = 0,
                                         cor_intercept_hu_slope = 0,
                                         cor_slope_hu_intercept = 0,
                                         cor_slope_hu_slope = 0,
                                         cor_hu_intercept_hu_slope = 0,
                                         sd_measures = NULL,
                                         sd_measures_hu = NULL,
                                         cor_measures = NULL,
                                         sigma_log = NULL,
                                         shape = NULL,
                                         dropout = 0,
                                         deterministic_dropout = FALSE,
                                         RR_cont = 1,
                                         OR_hu = 1,
                                         marginal = FALSE,
                                         EMA = FALSE,
                                         family = "lognormal", ...) {

    args <- as.list(match.call())[-1]

    args <- do.call(list, args)
    args$design <- NULL
    # dropout checks
    .check_dropout_arg(args$dropout)

    # defaults
    args$fixed_intercept_tx <- fixed_intercept_tx
    args$fixed_hu_intercept_tx <- fixed_hu_intercept_tx

    args$cor_intercept_slope <- cor_intercept_slope
    args$cor_intercept_hu_intercept <- cor_intercept_hu_intercept
    args$cor_intercept_hu_slope <- cor_intercept_hu_slope
    args$cor_slope_hu_intercept <- cor_slope_hu_intercept
    args$cor_slope_hu_slope <- cor_slope_hu_slope
    args$cor_hu_intercept_hu_slope <- cor_hu_intercept_hu_slope
    args$n3 <- 1 # 3-level not yet supported

    args$RR_cont <- RR_cont
    args$OR_hu <- OR_hu
    args$EMA <- EMA
    args$dropout <- dropout
    args$deterministic_dropout <- deterministic_dropout
    args$partially_nested <- FALSE
    save_call <- args
    save_call$design <- design

    pars <- expand.grid(args)
    pars$design <- "plcp_hurdle"

    ## Default T_end
    if(is.null(args$T_end)) pars$T_end <- pars$n1 - 1

    # Single or multi?
    pars <- .make_single_or_multi(pars, model = "hurdle")

    attr(pars, "call") <- save_call

    pars
}
