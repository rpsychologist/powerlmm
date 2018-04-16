# Design effect ----------------------------------------------------------------
approx_type1 <- function(bias) {
    pnorm(qnorm(0.025)* (1+bias)) * 2
}


#' Calculate the design effect and Type I errors
#'
#' This functions helps to evaluate the consequences of ignoring a random slope
#' at the cluster level.
#'
#' @param object A \code{plcp_3lvl}-object created by \code{\link{study_parameters}}
#'
#' @return A \code{data.frame} with the columns \code{n1, n2, n3, icc_slope,
#' var_ratio, DEFT}, and, \code{approx_type1}. The number of rows of the
#' \code{data.frame} will be equals to the number of
#' different combination of parameters values specified with \link{study_parameters}.
#' @seealso \code{\link{simulate.plcp}}
#' @export
#'
#' @details The design effect (DEFT) is the ratio of the standard error from the correct
#' three-level model to the standard error from the misspecified model omitting
#' the cluster-level random slope. The standard error for the misspecified model is
#' calculated by assuming that the cluster-level random slope variance is added
#' to the subject-level random slope.
#'
#' The approximate Type I error under the miss-specified model is also calculated.
#' The effect of wrongly ignoring a third-level random slope on the Type I errors, depends on
#' \code{n1, n2, n3, icc_slope}, and, \code{var_ratio}.
#'
#' @examples
#' paras <- study_parameters(n1 = 11,
#'                           n2 = 30,
#'                           n3 = 3,
#'                           T_end = 10,
#'                           icc_pre_subject = 0.5,
#'                           icc_pre_cluster = 0,
#'                           icc_slope = c(0.01,0.05, 0.1),
#'                           var_ratio = 0.02)
#'
#' get_DEFT(paras)
get_DEFT <- function(object) {
    UseMethod("get_DEFT")
}

#' @rdname get_DEFT
#' @export
get_DEFT.plcp_3lvl <- function(object) {
    if(any(object$partially_nested)) stop("DEFT is not yet implemented for partially nested designs. Use simulation, see '?simulate.plcp'")
    get_DEFT_3lvl(object)

}
get_DEFT_3lvl <- function(object) {
    UseMethod("get_DEFT_3lvl")
}
get_DEFT_3lvl.data.frame <- function(object) {
    res <- lapply(1:nrow(object), function(i) get_DEFT_3lvl.list(as.plcp(object[i,])))
    res <- do.call(rbind, res)
    res <- as.data.frame(res)
    res
}
get_DEFT_3lvl.list <- function(object) {
    p1 <- object
    p2 <- p1
    p2$sigma_subject_slope <- sqrt(p2$sigma_cluster_slope^2 +
                                       p2$sigma_subject_slope^2)
    p2$sigma_cluster_slope <- 0
    class(p2) <- gsub("3lvl", "2lvl", class(p2))


    x1 <- get_power(p1)
    x2 <- get_power(p2)

    se1 <- unlist(x1$se)
    se2 <- unlist(x2$se)

    DEFT <- se1/se2

    if(is.unequal_clusters(p1$n2)) warning("DEFT for designs with unequal clusters is highly experimental. Check your results using simulation, see ?simulate.plcp")

    data.frame(icc_slope = get_ICC_slope(p1),
               var_ratio = get_var_ratio(p1),
               DEFT = DEFT,
               approx_type1 = approx_type1((se2-se1)/se1))

}
get_rel_bias.paras <- function(n1, n2, n3, T_end, u1, v1, error, u0 = 0, u01 = 0, v0 = 0, v01 = 0) {

    sx <- Vectorize(var_T)(n1, T_end)
    n2 <- lapply(n2, function(x) mean(x))
    n2 <- unlist(n2)
    se1 <- sqrt(2*(error^2 + n1*u1^2*sx + n1*n2*v1^2*sx) / (n1*n2*n3*sx))
    se2 <- sqrt(2*(error^2 + n1*(u1^2+v1^2)*sx) / (n1*n2*n3*sx))

    # rel bias
    ICC_aov <- get_ICC_aov(u0 = u0,
                           u1 = u1,
                           u01 = u01,
                           t = T_end,
                           v0 = v0,
                           v1 = v1,
                           v01 = v01,
                           error=error)
    rel1 <- data.frame(n1 = n1,
                       n2 = n2,
                       n3 =  n3,
                       "approx_rel_bias" = (se2-se1)/se1,
                       "DEFT" = se1/se2,
                       "DEFT_aov" = sqrt(1 + (n2-1) * ICC_aov))

    data.frame(rel1,
               "ratio" = get_var_ratio(v1, u1, error),
               "ICC_slope" = get_ICC_slope(u1, v1),
               "ICC_aov" = ICC_aov, # icc at last T
               "sigma_subject_intercept" = u0,
               "approx_type1" = approx_type1(rel1$approx_rel_bias))

}

