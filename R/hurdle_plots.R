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
    x <- lapply(seq_along(args), function(i) {
        tmp <- args[[i]]
        tmp$var <- names(args)[i]

        tmp
    })
    x <- do.call(rbind, x)

    Q_long <- lapply(seq_along(args), function(i) {
        tmp <- reshape_eta_sum(args[[i]])
        tmp$var <- names(args)[i]

        tmp
    })
    Q_long <- do.call(rbind, Q_long)

    x$var <- factor(x$var, levels = c("overall", "positive", "hurdle"))

    # Get limits
    if(RE) {
        lims <- lapply(names(args), function(x) {
            tmp <- Q_long[Q_long$var == x, ]
            data.frame(var = x,
                       mean = c(min(tmp$min), max(tmp$max)),
                       treatment = "Treatment",
                       time = 0)
        })
    } else {
        lims <- lapply(names(args), function(y) {
            tmp <- x[x$var == y, ]
            data.frame(var = y,
                       mean = c(min(tmp$mean, tmp$Q50), max(tmp$mean, tmp$Q50)),
                       treatment = "Treatment",
                       time = 0)
        })
    }
    lims <- do.call(rbind, lims)

    # use same limits for 'overall' and 'positive'
    if(all(lims$var %in% c("overall", "positive"))) {
        tmp <- lims[lims$var %in% c("overall", "positive"), ]
        lims[lims$var == "overall", "mean"] <- c(min(tmp$mean), max(tmp$mean))
        lims[lims$var == "positive", "mean"] <- c(min(tmp$mean), max(tmp$mean))
    }

    lims$var <- factor(lims$var, levels = c("overall", "positive", "hurdle"))

    if(type == "dropout") {
        .plot_dropout(object$paras)
    } else if(type %in% c("trend", "trend_dropout")) {
        if(type == "trend_dropout") {
            check_installed("gridExtra")
            pd <- .plot_dropout(object$paras)
        }

        p2 <- .plot_marg(x = x,
                         Q_long = Q_long,
                         RE = RE,
                         ymin = ymin,
                         ymax = ymax) +
            geom_blank(data = lims) +
            labs(linetype = "", color = "", title = "Subject level") +
            facet_wrap(~var, scales = "free")
        if(type == "trend_dropout") {
            gridExtra::grid.arrange(p2, pd, ncol=1)
            return(invisible(list("trend" =  list("subject" = p2),
                                  "dropout" = pd)))
        } else {
            plot(p2)
            return(invisible(p2))
        }

    } else if(type %in% c("post_diff", "post_ratio", "post_ratio_diff")) {
        .plot_diff_marg(object, type = type, ...)
    }

}
