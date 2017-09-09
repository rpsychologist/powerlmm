
make_dropout_indicator_deterministic <- function(time, retention, tot_n) {
    d <- data.frame(id = rep(1:tot_n, each = length(time)),
                    time = rep(time, tot_n),
                    missing = 0)

    ids <- unique(d$id)

    for(i in seq_along(retention)) {
        n_drop <- round((length(ids)) - (retention[i] * tot_n))
        dropouts <- sample(ids, n_drop)
        ids <- ids[which(!ids %in% dropouts)]
        d[d$id %in% dropouts & d$time >= time[i], "missing"] <- 1
    }

    d
}


dropout_per_time <- function(time, retention, tot_n) {
    id <- 1:tot_n
    missing <- rep(0, each=retention * tot_n)
    missing <- c(missing, rep(1, tot_n - length(missing)))
    data.frame(id, time, missing = missing)
}

# stochastic dropout
make_dropout_indicator_stochastic <- function(time, retention, tot_n) {

    ps <-  abs(diff(retention, 1))
    ps <- c(ps, tail(retention, 1))
    miss_ind <- rmultinom(tot_n, 1, prob = ps)

    d <- data.frame(id = rep(1:tot_n, each = length(time)),
                    time = rep(time, tot_n),
                    miss_ind = c(miss_ind),
                    missing = 1)

    for(i in unique(d$id)) {
        x <- d[d$id == i, ]
        missing <- x$missing
        missing[1:which(x$miss_ind == 1)] <- 0
        d[d$id == i, "missing"] <- missing
    }

    d
}


dropout_process <- function(dropout, paras) {

    if(paras$deterministic_dropout) {
        make_dropout_indicator <- make_dropout_indicator_deterministic
    } else {
        make_dropout_indicator <- make_dropout_indicator_stochastic
    }

    time <- get_time_vector(paras)
    if(is.function(dropout[[1]])) {
        retention <- dropout[[1]](time)
    } else retention <- 1 - dropout
    n2 <- unlist(paras$n2)
    if(length(n2) > 1) {
        tot_n <- sum(n2)
    } else {
        tot_n <- n2 * paras$n3
    }

    if(is.list(retention)) {
        d_a <- make_dropout_indicator(time, retention[[1]], tot_n)
        d_a$treatment <- "A"

        d_b <- make_dropout_indicator(time, retention[[2]], tot_n)
        d_b$treatment <- "B"

        res <- list(d_a, d_b)
    } else {
        d <- make_dropout_indicator(time, retention, tot_n)

        res <- d
    }

    res
}


#' Get the amount of dropout
#'
#' @param object An object created by \code{\link{study_parameters}}
#' @param n The \emph{n}-th dataset to use for objects with multiple designs.
#' @param ... Optional arguments.
#'
#' @return A \code{data.frame} with the proportion of dropout per time point
#' and treatment condition.
#' @seealso \code{\link{dropout_manual}}, \code{\link{dropout_weibull}}
#' @export
#'
#' @examples
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout_weibull(proportion = 0.3, rate = 3),
#'                       cohend = -0.8)
#' get_dropout(p)

#' @export
get_dropout <- function(object, ...) {
    UseMethod("get_dropout")
}

#' @export
get_dropout.plcp <- function(object, ...) {
    paras <- object
    dropout <- paras$dropout
    time <- get_time_vector(paras)
    dropout_cc <- dropout
    dropout_tx <- dropout

    if(is.per_treatment(paras$dropout)) {
        dropout_cc <- unlist(dropout[[1]]$control)
        dropout_tx <- unlist(dropout[[1]]$treatment)
    } else {
        dropout_cc <- unlist(dropout)
        dropout_tx <- unlist(dropout)
    }

    if(is.function(dropout_cc[[1]])) {
        dropout_cc <- 1-dropout_cc[[1]](time)
    }
    if(is.function(dropout_tx[[1]])) {
        dropout_tx <- 1-dropout_tx[[1]](time)
    }

    data.frame(time = time,
               control = dropout_cc,
               treatment = dropout_tx)
}
#' @rdname get_dropout
get_dropout.plcp_multi <- function(object, n = 1, ...) {
    get_dropout.plcp(object[n,])
}
## weibull dropout

#' Use the Weibull distribution to specify the dropout process
#'
#' Used as input to the \code{dropout}-argument in \code{\link{study_parameters}}
#'
#' @param proportion Total proportion of subjects that have dropped out
#' at the last time point. Must be less than 1.
#' @param rate Indicates the "shape" of the dropout process, if > 1 then dropout is
#' concentrated at the end of the study, if \code{rate} < 1 more dropout occurs at
#' the beginning of the study. If \code{rate} == 1 the risk of dropout is constant.
#'
#' @details N.B a constant (rate = 1) hazard of dropout does not mean dropout
#' is linear over time. It means that the risk of dropping out at the next
#' time point is constant over the study period.
#'
#' @seealso \code{\link{dropout_manual}}, \code{\link{per_treatment}}
#'
#' @references Galbraith, S., Stat, M., & Marschner, I. C. (2002).
#' Guidelines for the design of clinical trials with longitudinal outcomes.
#' \emph{Controlled clinical trials, 23}(3), 257-273.
#'
#' @return A \code{plcp_weibull} named \code{list}, with the first element containing the
#' dropout \code{function}.
#' @export
#'
#' @examples
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout_weibull(proportion = 0.3, rate = 3),
#'                       cohend = -0.8)
#'
#' get_dropout(p)
#' plot(p, plot = 2)
#'
#' # Different per treatment
#' tx <- dropout_weibull(proportion = 0.3, rate = 3)
#' cc <- dropout_weibull(proportion = 0.3, rate = 1/3)
#' dropout <- per_treatment(control = cc,
#'                          treatment = tx)
#'
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout,
#'                       cohend = -0.8)
#'
#' plot(p, plot = 2)
#'
#' # Compare power for different dropout amounts
#' dropout <- c(dropout_weibull(proportion = 0.3, rate = 3),
#'              dropout_weibull(proportion = 0.5, rate = 3),
#'              dropout_weibull(proportion = 0.5, rate = 1/3))
#'
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout,
#'                       cohend = -0.8)
#'
#' get_power(p)
dropout_weibull <- function(proportion, rate) {
    if(proportion == 1) stop("'proportion' must be less than 1")
    if(proportion < 0) stop("'proportion' can't be negative")
    if(rate <= 0) stop("'rate' must be greater than zero")
    f <- function(time) {
        time <- time/max(time)
        (1 -proportion)^(time^rate)
    }

    f <- list("dropout_weibull" = f)
    class(f) <- append(class(f), "plcp_weibull")

    f
}

#' Manually specify dropout per time point
#'
#' Used as input to the \code{dropout}-argument in \code{\link{study_parameters}}.
#'
#' @param ... The proportion of dropout per time point, either as
#' a vector of length \code{n1}, or \code{n1} individual numeric arguments,
#' see \emph{Details}.
#' @details Specifying dropout manually requires that the dropout
#' is 0 at the first time point. Moreover, dropout can't decrease over time and
#' can never be 1.
#'
#' @return A list of class \code{plcp_dropout_manual}
#' @seealso \code{\link{dropout_weibull}}, \code{\link{per_treatment}}
#' @export
#'
#' @examples
#' dropout <- dropout_manual(0, 0, 0, 0, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.45)
#'
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout,
#'                       cohend = -0.8)
#' plot(p, plot = 2)
#' get_power(p)
#'
#' # Can also use a vector as input
#' dropout <- dropout_manual(seq(0, 0.5, length.out = 11))
#' p <- study_parameters(n1 = 11,
#'                       n2 = 5,
#'                       n3 = 6,
#'                       T_end = 10,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0,
#'                       var_ratio = 0.03,
#'                       icc_slope = 0.05,
#'                       dropout = dropout,
#'                       cohend = -0.8)
#' plot(p, plot = 2)
#' get_power(p)
#'
#' \dontrun{
#' # Decreasing dropout will throw an error
#' dropout_manual(0, 0.1, 0.1, 0.2, 0.1)
#'
#' # Dropout at the first time point will throw an error
#' dropout_manual(0.1, 0.1, 0.1, 0.2, 0.2)
#' }

dropout_manual <- function(...) {
    x <- list(dropout_manual = ...)

    tmp <- unlist(x)
    if(any(tmp >= 1)) stop("Values must be less than 1")
    if(any(tmp < 0)) stop("Values can't be negative")
    if(tmp[1] != 0) stop("Dropout must be 0 at the first time point")
    if(any(sign(diff(tmp)) == -1)) stop("Dropout can't be decreasing over time")


    class(x) <- "plcp_dropout_manual"
    x <- list(dropout_manual = x)
    #class(x) <- "unequal_clusters"
    x
}




