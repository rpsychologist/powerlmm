plot.plcp_marginal_hurdle <- function(object,
                                      type = "trend",
                                      outcome = c("overall", "positive", "hurdle"),
                                      RE = TRUE, RE_level = 2, hu = FALSE, trim = c(0,1),
                                      ...) {
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

        if(any(RE_level %in% c(1,3))) message("'RE_level' is ignored, level 2 is shown.")
        # trend <- .make_nested_trend(object = object,
        #                             RE = RE,
        #                             RE_level = 2,
        #                             var1 = "y",
        #                             var2 = "y_overall",
        #                             level1_func = .sample_level1_nested_hurdle,
        #                             ...)

        # trend <- .make_nested_trend(object = object,
        #                             RE = RE,
        #                             RE_level = RE_level,
        #                             var1 = "y",
        #                             var2 = "y_overall",
        #                             level1_func = .sample_level1_nested_hurdle,
        #                             ...)
        #
        # trend$x$var <- "Subject_overall"
        # trend$Q_long$var <- "Subject_overall"
        #
        #
        # trend1 <- .make_nested_trend(object = object,
        #                             RE = RE,
        #                             RE_level = RE_level,
        #                             var1 = "y",
        #                             var2 = "y_positive",
        #                             level1_func = .sample_level1_nested_hurdle,
        #                             ...)
        # trend1$x$var <- "Subject_positive"
        # trend1$Q_long$var <- "Subject_positive"
        #
        # trend2 <- .make_nested_trend(object = object,
        #                              RE = RE,
        #                              RE_level = RE_level,
        #                              var1 = "y",
        #                              var2 = "hu_prob",
        #                              level1_func = .sample_level1_nested_hurdle,
        #                              ...)
        #
        # trend2$x$var <- "Hu prob"
        # trend2$Q_long$var <- "Hu prob"
        #
        # trend$x <- rbind(trend$x, trend1$x, trend2$x)
        # trend$Q_long <- rbind(trend$Q_long, trend1$Q_long, trend2$Q_long)

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
                                     var1 = "y",
                                     var2 = "y_overall",
                                     level1_func = .sample_level1_nested_hurdle,
                                     ...)

         #res <- subset(res, y > trend$lims$mean[1] & y < trend$lims$mean[2])


         # TODO: trim per facet, usefull to see range of data
         res <- lapply(unique(res$var), function(x) {
             tmp <- res[res$var == x, ]
             y <-  quantile(tmp$y, trim)
             lims <- data.frame(y = y, var = x, treatment = "Treatment", time = 0)

             subset(tmp, y >= min(lims$y) & y <= max(lims$y))

         })
         res <- do.call(rbind, res)

         res_dens <- res %>%
             filter(y > 0) %>%
             group_by(time, treatment, var) %>%
             do(ggplot2:::compute_density(.$y, NULL)) %>%
             filter(x > 0)

        #res <- subset(res, y >= min(lims$y) & y <= max(lims$y))

         max_count <- res_dens %>%
             filter(var == "Within-subject") %>%
             summarise(count = max(count))

         res0 <- res %>%
             filter(var == "Within-subject") %>%
             group_by(time, treatment, var) %>%
             summarise(p = mean(near(y, 0)))

         step <- mean(diff(unique(res0$time)))
         res0 <- res0 %>%
             group_by(treatment) %>%
             mutate(xend = time + step * p)

         sub_clust <- subset(trend$x, var %in% c("Subject", "Cluster"))

         p <- ggplot(res, aes(x = y,
                         y = time,
                         group = interaction(time, treatment, var),
                         fill = treatment,
                         color = treatment)) +
             # ggridges::geom_density_ridges(data = res_dens,
             #                               scale = 1,
             #                               stat = "identity",
             #                               aes(x = x, height = count, color = NULL),
             #                               color = NA,
             #                               alpha = 0.75,
             #                               trim = TRUE) +
             ggridges::geom_density_ridges(data = filter(res, y > 0),
                                           scale = 1,
                                           stat = "density",
                                           aes(height = ..count.., color = NULL),
                                           color = NA,
                                           alpha = 0.75,
                                           draw_baseline = FALSE,
                                           trim = TRUE) +
             #ggridges::stat_density_ridges(data = res, quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
             geom_segment(data = res0,
                          aes(y = time, yend = xend, x = 0, xend = 0, color = treatment),
                          alpha = 0.5,
                          stat = "identity", size = 2) +

             # ggridges::geom_density_ridges(data = subset(res, y == 0),
             #                               stat = "binline",
             #                               scale = 0.95,
             #                               binwidth = 1,
             #                               aes(height = ..count.., color = NULL),
             #                               color = alpha("white", 0.5),
             #                               alpha = 0.75,
             #                               size = 0.5,
             #                               draw_baseline = FALSE) +
             geom_path(data = trend$x,
                       aes(x = mean,
                           y = time,
                           linetype = "mean",
                           fill = NULL,
                           group = interaction(treatment, var)),
                       size = 1)
             #TODO
             # geom_segment(data = trend$x,
             #              aes(x = mean,
             #                  xend = mean,
             #                  y = time,
             #                  yend = time+1,
             #                  linetype = "mean",
             #                  fill = NULL,
             #                  group = interaction(treatment, time, var)),
             #              size = 1) +

            if(nrow(sub_clust) > 0) {
                p <- p +
                    geom_path(data = sub_clust,
                              aes(x = Q50,
                                  y = time,
                                  linetype = "median",
                                  fill = NULL,
                                  group = interaction(treatment, var)), size = 1)
            }

             #geom_blank(data = trend$lims, aes(x = mean, y = time)) +
             p +
             coord_flip() +
             theme_minimal() +
             facet_wrap(~var, ncol = 2, scales = NULL) +
             scale_fill_manual(values = c("#30394F", "#6ACEEB")) +
             scale_color_manual(values = c("#222f3e", "#2980b9"))


    } else if(type %in% c("post_diff", "post_ratio", "post_ratio_diff")) {
        .plot_diff_marg(object, type = type, ...)
    }

}




# Helpers -----------------------------------------------------------------
## ggplot2:::compute_density
compute_density <- function (x, w = NULL, from, to, bw = "nrd0", adjust = 1, kernel = "gaussian",
          n = 512) {
    nx <- length(x)
    if (is.null(w)) {
        w <- rep(1/nx, nx)
    }
    if (nx < 2) {
        warning("Groups with fewer than two data points have been dropped.",
                call. = FALSE)
        return(data.frame(x = NA_real_, density = NA_real_, scaled = NA_real_,
                          count = NA_real_, n = NA_integer_))
    }
    dens <- stats::density(x, weights = w, bw = bw, adjust = adjust,
                           kernel = kernel, n = n, from = from, to = to)
    data.frame(x = dens$x,
               density = dens$y,
               count = dens$y * nx)
}



