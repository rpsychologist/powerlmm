
eta_sum <- function(x) {
    # return summaries of marginal ests distribution
    cbind("mean" = mean(x),
          "sd" = sd(x),
          "Q0.5" = quantile(x, probs = 0.005),
          "Q2.5" = quantile(x, probs = 0.025),
          "Q10" = quantile(x, probs = 0.10),
          "Q25" = quantile(x, probs = 0.25),
          "Q50" = median(x),
          "Q75" = quantile(x, probs = 0.75),
          "Q90" = quantile(x, probs = 0.90),
          "Q97.5" = quantile(x, probs = 0.975),
          "Q99.5" = quantile(x, probs = 0.995)

    )

}


trans_eta <- function(x, var, d) {
    # helper to display marginal ests
    out <- do.call(rbind, x[, var])
    out <- as.data.frame(out)
    out <- cbind(data.frame(var = var,
                            treatment = d$treatment,
                            time = d$time),
                 out)
    rownames(out) <- NULL
    out
}


trans_post_ps <- function(x, hu = FALSE) {
    # transform percentile ES to long foramt
    post_ps <- x[vapply(x, is.data.frame, logical(1))]

    ind <- vapply(post_ps, function(d) all(d$treatment == 0), logical(1))
    post_ps_c <- post_ps[[which(ind)]]
    post_ps_tx <- post_ps[[which(!ind)]]

    if(hu) {
        post_ps <- data.frame(percentile = post_ps_c$percentile,
                              diff = post_ps_tx$value - post_ps_c$value,
                              OR = get_OR(post_ps_tx$value, post_ps_c$value)
        )
    } else {
        post_ps <- data.frame(percentile = post_ps_c$percentile,
                              diff = post_ps_tx$value - post_ps_c$value,
                              ratio = post_ps_tx$value / post_ps_c$value)
    }


    list(control = post_ps_c,
         treatment = post_ps_tx,
         effect = post_ps)
}
