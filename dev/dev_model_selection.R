p <- study_parameters(n1 = 11,
                      n2 = 20,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      partially_nested = TRUE,
                      var_ratio = 0.03)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)


res <- simulate(p, formula = f, nsim = 1000, satterthwaite = TRUE, cores = 10, CI = FALSE)
summary(res)

winners <- step.plcp_sim(res, alpha = 0.5)


# TODO: write function to use 'winners' to pick parameters from winning models and summarise

nsim <- res$nsim
x <- list("RE" = vector("list", nsim),
          "FE" = vector("list", nsim),
          "tot_n" = vector("numeric", nsim),
          "convergence" = vector("numeric", nsim))

get_sim_para <- function(i, effect, res, mod) {
    x <- res$res[[mod]][[effect]]
    x <- x[x$sim == i, ]
    x$mod <- mod
    x
}

for(i in seq_along(winners)) {
    mod <- winners[i]
    x$RE[[i]] <- get_sim_para(i, "RE", res, mod)
    x$FE[[i]] <- get_sim_para(i, "FE", res, mod)
    x$tot_n[i] <- res$res[[mod]][["tot_n"]][i]
    x$convergence[i] <- res$res[[mod]][["convergence"]][i]
}
x$RE <- do.call(rbind, x$RE)
x$FE <- do.call(rbind, x$FE)

res2 <- res
res2$res$test <- x
summary(res2)

## TODO: enable to summarise a parameter as a function of step.plcp_sum(alpha = x).
## useful to plot e.g. type I errors for time:treatment as a function of LRT alpha level.

res2$res$test$FE %>%
    group_by(parameter) %>%
    summarise(mean(estimate),
              mean(se),
              sd(estimate))

