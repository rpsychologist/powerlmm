
#' @param p1 prob tx
#' @param p0 prob control
get_OR <- function(p1, p0) {
    (p1/(1-p1)) / (p0 / (1-p0))
}

get_ICC_pre_subjects.plcp_hurdle <- function(object, ...) {
    numeric(0)
}
get_ICC_pre_subjects.plcp_multi_hurdle <- function(object, ...) {
    numeric(nrow(object))
}
get_ICC_pre_clusters.plcp_hurdle <- function(object, ...) {
    numeric(0)
}
get_ICC_pre_clusters.plcp_multi_hurdle <- function(object, ...) {
    numeric(nrow(object))
}

get_ICC_slope.plcp_hurdle  <- function(object, ...) {
    numeric(0)
}
get_ICC_slope.plcp_multi_hurdle <- function(object, ...) {
    numeric(nrow(object))
}
get_var_ratio.plcp_hurdle  <- function(object, ...) {

    warning("Not yet implemented for 'hurdle' models")
    numeric(0)
}
get_var_ratio.plcp_multi_hurdle <- function(object, ...) {
    numeric(nrow(object))
}


# CTP ---------------------------------------------------------------------

get_overall_hurdle <- function(p) {

    if(p$marginal) {
        if(p$family == "gamma") {
            cc_mu <- p$fixed_intercept + p$fixed_slope * p$T_end
            tx_mu <- p$fixed_intercept + (p$fixed_slope) * p$T_end + log(p$RR_cont)

            tx <- exp(tx_mu)
            cc <- exp(cc_mu)

        } else if(p$family == "lognormal") {
            stop("not implemented")
        }
    } else {
        if(p$family == "gamma") {
            cc_logit <- p$fixed_hu_intercept + p$fixed_hu_slope * p$T_end
            tx_logit <- p$fixed_hu_intercept + (p$fixed_hu_slope) * p$T_end + log(p$OR_hu)
            cc_mu <- p$fixed_intercept + p$fixed_slope * p$T_end
            tx_mu <- p$fixed_intercept + (p$fixed_slope) * p$T_end + log(p$RR_cont)

            tx <- plogis(-tx_logit) * exp(tx_mu)
            cc <- plogis(-cc_logit) * exp(cc_mu)

        } else if(p$family == "lognormal") {
            stop("not implemented")
        }
    }


    list("tx" = tx,
         "cc" = cc,
        "diff" = tx - cc,
        "ratio" = tx/cc)
}


# Plot --------------------------------------------------------------------

plot_hurdle_time <- function(x) {

    tmp <- x[, !colnames(x) %in% c("var","mean", "sd", "Q50")]
    tmp <- lapply(list(c("Q0.5", "Q99.5"),
                       c("Q2.5", "Q97.5"),
                       c("Q10", "Q90"),
                       c("Q25", "Q75")),
                  function(x) {
                      d <- tmp[, c("treatment","time", x)]
                      data.frame(treatment = d$treatment,
                                 time = d$time,
                                 min = d[, x[1]],
                                 max = d[, x[2]],
                                 width = switch(x[1],
                                                "Q0.5" = 0.99,
                                                "Q2.5" = 0.95,
                                                "Q10" = 0.8,
                                                "Q25" = 0.5)
                      )
                  })
    tmp <- do.call(rbind, tmp)
    tmp$width <- factor(tmp$width,
                        levels = rev(c(0.5, 0.8, 0.95, 0.99)))

    ggplot(x, aes(time, mean, group = treatment)) +
        geom_ribbon(data = tmp, aes(ymin = min, ymax = max, y = NULL, x = time, group = width, fill = width), alpha = 0.75) +
        geom_line(aes(color = "mean", linetype = "mean", fill = NULL), size = 1) +
        geom_line(aes(y = Q50, color = "median", linetype = "median", fill = NULL), size = 1) +
        geom_point(aes(y = Q50), color = "red") +
        scale_color_manual(values = c("median" = "red", "mean" = "red")) +
        scale_linetype_manual(values = c("median" = "solid", "mean" = "dotted")) +
        labs(linetype = "", color = "") +
        guides(color = guide_legend(override.aes = list(fill = NA))) +
        facet_wrap(~treatment, ncol = 2) +
        scale_fill_brewer() +
        theme_minimal()
}

plot_hurdle_diff <- function(x, hu = FALSE, fixed_overall = NULL) {

    if(hu) {
        tmp <- x$post_hu_ps$effect
        ES <- x$post_hu[x$post_hu$var == "marg_hu_post_diff", "est"]
        ES_med <- x$post_hu[x$post_hu$var == "median_hu_post_diff", "est"]
        ES_ratio <- x$post_hu[x$post_hu$var == "marg_OR", "est"]
        ES_ratio_med <- x$post_hu[x$post_hu$var == "median_OR", "est"]
    } else {
        tmp <- x$post_ps$effect
        ES <- x$post[x$post$var == "marg_post_diff", "est"]
        ES_med <- x$post[x$post$var == "median_post_diff", "est"]
        ES_ratio <- x$post[x$post$var == "marg_RR", "est"]
        ES_ratio_med <- x$post[x$post$var == "median_RR", "est"]
    }

    tmp$fill <- ifelse(tmp$percentile == 0.5, "median", "other")
    tmp$fill[which.min(abs(tmp$diff - ES))] <- "mean"


    p0 <-  ggplot(tmp, aes(percentile, diff, fill = fill)) +
        geom_histogram(stat = "identity", color = "white", fill = "#3498db", alpha = .75) +
        #geom_hline(yintercept = 0, linetype = "solid", size = 0.75) +
        geom_hline(yintercept = ES, linetype = "dotted", alpha = 0.75, size = 0.75) +
        geom_hline(yintercept = ES_med, linetype = "dashed", alpha = 0.75, size = 0.75) +
        scale_y_continuous(sec.axis = sec_axis(~ ., breaks = c(ES_med, ES),
                                               labels = c(paste(round(ES_med, 2), " (median)"),
                                                          paste(round(ES, 2), " (mean)")
                                               )
        )
        ) + theme_minimal() +
        theme(legend.position = "none")


    # Ratio
    if(hu) tmp$ratio <- tmp$OR
    if(abs(ES_ratio_med - ES_ratio) < .Machine$double.eps^0.5) {
        breaks <- ES_ratio
        labels <- paste(round(ES_ratio, 2), " \n(mean,\nmedian)")
    } else {
        breaks <- c(ES_ratio_med, ES_ratio)
        labels <- c(paste(round(ES_ratio_med, 2), " (median)"),
                    paste(round(ES_ratio, 2), " (mean)")
                    )
    }

    p1 <- ggplot(tmp, aes(percentile, ratio)) +
        geom_histogram(stat = "identity", color = "white", fill = "#3498db", alpha = .75) +
        #geom_hline(yintercept = 0, linetype = "dotted", size = 0.75) +
        geom_hline(yintercept = ES_ratio, linetype = "dotted", size = 0.75) +
        geom_hline(yintercept = ES_ratio_med, linetype = "dashed", size = 0.75) +
        scale_y_continuous(sec.axis = sec_axis(~ ., breaks = breaks,
                                               labels = labels
                                               )) +
        theme_minimal() +
        theme(legend.position = "none")

    if(!is.null(fixed_overall)) {
        p0 <- p0 + geom_hline(yintercept = fixed_overall$diff, color = "#e74c3c")
        p1 <- p1 + geom_hline(yintercept = fixed_overall$ratio, color = "#e74c3c")
    }

    gridExtra::grid.arrange(p0, p1)

}

plot_hurdle_probs <- function(x) {
    # cc
    tmp <- x$post_hu_ps$control
    ES <- x$post_hu[x$post_hu$var == "marg_hu_post_cc", "est"]
    tmp$fill <- ifelse(tmp$percentile == 0.5, "median", "other")
    tmp$fill[which.min(abs(tmp$value - ES))] <- "mean"

    # tx
    tmp_tx <- x$post_hu_ps$treatment
    ES_tx <- x$post_hu[x$post_hu$var == "marg_hu_post_tx", "est"]
    tmp_tx$fill <- ifelse(tmp_tx$percentile == 0.5, "median", "other")
    tmp_tx$fill[which.min(abs(tmp_tx$value - ES_tx))] <- "mean"

    tmp <- rbind(tmp, tmp_tx)

    c_p_mean <- x$post_hu[x$post_hu$var == "marg_hu_post_cc", "est"]
    c_p_median <-  x$post_hu[x$post_hu$var == "median_hu_post_cc", "est"]

    tx_p_mean <- x$post_hu[x$post_hu$var == "marg_hu_post_tx", "est"]
    tx_p_median <-  x$post_hu[x$post_hu$var == "median_hu_post_tx", "est"]

    median_x_pos <- mean(tmp[tmp$fill == "median", "percentile"])
    mean_x_pos <- mean(tmp[tmp$fill == "mean", "percentile"])

    median_y_lab_pos <- mean(tmp[tmp$fill == "median", "value"])
    mean_y_lab_pos <- mean(tmp[tmp$fill == "mean", "value"])

    ggplot(tmp, aes(percentile, value, group = treatment, color = fill)) +
        geom_point() +
        geom_hline(yintercept = tx_p_mean, color = "cyan", alpha = 0.5) +
        geom_hline(yintercept = c_p_mean, color = "cyan", alpha = 0.5) +

        geom_hline(yintercept = tx_p_median, color = "red", alpha = 0.5) +
        geom_hline(yintercept = c_p_median, color = "red", alpha = 0.5) +
        annotate("segment",
                 x = mean_x_pos,
                 xend = mean_x_pos,
                 y = tx_p_mean,
                 yend = c_p_mean) +
        annotate("segment",
                 x = median_x_pos,
                 xend = median_x_pos,
                 y = tx_p_median,
                 yend = c_p_median) +
        annotate("label", x = mean_x_pos, y = mean_y_lab_pos, label = paste0("Mean OR: ", round(get_OR(tx_p_mean, c_p_mean), 2)), hjust = -1) +
        annotate("label", x = median_x_pos, y = median_y_lab_pos, label = paste0("Median OR: ", round(get_OR(tx_p_median, c_p_median), 2)), hjust = -1) +
        labs(y = "Pr(Y = 0)") +
        scale_color_manual(values = c("median" = "red", "mean" = "cyan", "other" = "gray")) +
        theme_minimal() +
        theme(legend.position = "none")
}
