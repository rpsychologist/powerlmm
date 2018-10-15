plot.plcp_marginal_hurdle <- function(object, type = "trend", outcome = c("overall", "positive", "hurdle"),
                                      RE = TRUE, RE_level = 2, hu = FALSE, ...) {
    check_installed("ggplot2")
    # lvl 2
    p_overall <- NULL
    p_positive <- NULL
    p_hu <- NULL
    if(any(outcome == "overall")) {
        p_overall <- object$y_overall
    }
    if(any(outcome == "positive")) {
        p_positive <- object$y_positive
    }
    if(any(outcome == "hurdle")) {
        p_hu <- object$hu_prob
    }

    args <- list("overall" = p_overall,
                 "positive" = p_positive,
                 "hurdle" = p_hu)
    args <- args[!vapply(args, is.null, logical(1))]
    x <- .rbind_lists(args)
    Q_long <- .rbind_lists(args, func = reshape_eta_sum)

    x$var <- factor(x$var, levels = c("overall", "positive", "hurdle"))
    x$color <- "Subject"
    Q_long$color <- "Subject"

    # Get limits
    if(RE) {
        lims <- .get_facet_lims(d = Q_long,
                                var_names = names(args),
                                min_cols = "min",
                                max_cols = "max")
    } else {
        lims <- .get_facet_lims(d = x,
                                var_names = names(args),
                                min_cols = c("mean", "Q50"),
                                max_cols = c("mean", "Q50"))
    }

    # use same limits for 'overall' and 'positive'
    if(all(c("overall", "positive") %in% lims$var)) {
        tmp <- lims[lims$var %in% c("overall", "positive"), ]
        lims[lims$var == "overall", "mean"] <- c(min(tmp$mean), max(tmp$mean))
        lims[lims$var == "positive", "mean"] <- c(min(tmp$mean), max(tmp$mean))
    }

    lims$var <- factor(lims$var, levels = c("overall", "positive", "hurdle"))

    if(type == "dropout") {
        .plot_dropout(object$paras)
    } else if(type %in% "trend") {

        .plot_marg(x = x,
                   Q_long = Q_long,
                   RE = RE,
                   ymin = ymin,
                   ymax = ymax) +
            geom_blank(data = lims) +
            labs(linetype = "", color = "", title = "Subject level") +
            facet_wrap(treatment~var, scales = "free")
    } else if(type == "trend_ridges") {

         res <- .mu_vec_to_long(object,
                                RE_level = RE_level,
                                var1 = "mu1_vec",
                                var2 = "mu_overall_vec",
                                tx_var = "y_overall",
                                level1_func = .sample_level1_nested_hurdle,
                                ...)
         trend <- .make_nested_trend(object = object,
                                     RE = RE,
                                     RE_level = RE_level,
                                     var1 = "y_overall",
                                     var2 = "y_overall",
                                     level1_func = .sample_level1_nested_hurdle,
                                     ...)

         #res <- subset(res, y > trend$lims$mean[1] & y < trend$lims$mean[2])


        ggplot(res, aes(x = y,
                        y = time,
                        group = interaction(time, treatment, var),
                        fill = treatment, color = treatment)) +
            ggridges::geom_density_ridges(scale = 1,
                                          stat = "density",
                                          binwidth = 1,
                                          aes(height = ..count..),
                                          #binwidth = 1,
                                          #rel_min_height = 0.01,
                                          #color = alpha("black", 0.5),
                                          alpha = 0.33,
                                          size = 0.3,
                                          trim = FALSE) +
            geom_line(data = trend$x,
                      aes(x = mean,
                          y = time,
                          linetype = "mean",
                          fill = NULL,
                          group = interaction(treatment, var)),
                      size = 1) +
            geom_line(data = trend$x,
                      aes(x = Q50,
                          y = time,
                          linetype = "median",
                          fill = NULL,
                          group = interaction(treatment, var)),
                      size = 1) +
            coord_flip() +
            theme_minimal() +
            facet_wrap(~var, ncol = 2)

    } else if(type %in% c("post_diff", "post_ratio", "post_ratio_diff")) {
        .plot_diff_marg(object, type = type, ...)
    }

}
