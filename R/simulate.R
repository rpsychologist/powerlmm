#' Generate a data set using a \code{study_parameters}-object
#'
#' @param paras An object created by \code{\link{study_parameters}}
#' @param n Optional; specifies which row \code{n} should be used if \code{object}
#'  is a \code{data.frame} containing multiple setups.
#'
#' @return A \code{data.frame} with the simulated data in long form. Column
#' \code{y} is the outcome, \code{time} the time variable, \code{subject} is
#' the subject-level id variable, from 1, ..., the total number of subjects. \code{cluster} is
#' the cluster-level id variable, from 1, ..., to the total number of clusters.
#' @export
#'
#' @examples
#' p <- study_parameters(n1 = 11,
#'                       n2 = 10,
#'                       n3 = 4,
#'                       T_end = 10,
#'                       fixed_intercept = 37,
#'                       fixed_slope = -0.65,
#'                       sigma_subject_intercept = 2.89,
#'                       sigma_cluster_intercept = 0.6,
#'                       icc_slope = 0.1,
#'                       var_ratio = 0.03,
#'                       sigma_error = 1.5,
#'                       cor_subject = -0.5,
#'                       cor_cluster = 0,
#'                       cohend = 0.5)
#'
#' d <- simulate_data(p)
simulate_data <- function(paras, n = 1) {
    UseMethod("simulate_data")
}

#' @rdname simulate_data
#' @export
simulate_data.plcp <- function(paras, n = NULL) {
    if (is.data.frame(paras))
        paras <- as.list(paras)
    slope_diff <- get_slope_diff(paras) / paras$T_end
    paras$cohend <- NULL
    tmp <- prepare_paras(paras)
    paras <- tmp$control
    paras_tx <- tmp$treatment

    paras_tx$fixed_slope <- paras_tx$fixed_slope + slope_diff

    paras_tx$partially_nested <- NULL
    paras$partially_nested <- NULL
    paras$allocation_ratio <- NULL
    paras_tx$allocation_ratio <- NULL


    d_tx <- simulate_3lvl_data(paras_tx)
    d_c <- simulate_3lvl_data(paras)

    # drop outs
    if (is.list(paras$dropout) |
        is.function(paras$dropout) | length(paras$dropout) > 1) {
        miss_tx <- dropout_process(unlist(paras_tx$dropout), paras_tx)
        miss_tx <- miss_tx$missing

        miss_c <- dropout_process(unlist(paras$dropout), paras)
        miss_c <- miss_c$missing

        d_tx$miss <- miss_tx
        d_tx$y <- ifelse(d_tx$miss == 1, NA, d_tx$y)

        d_c$miss <- miss_c
        d_c$y <- ifelse(d_c$miss == 1, NA, d_c$y)
    }



    # combine
    d_tx$treatment <- 1
    d_c$treatment <- 0
    d_c$subject <- d_c$subject + max(d_tx$subject)
    d_tx$cluster <- d_tx$cluster + max(d_c$cluster)
    dt <- rbind(d_c, d_tx)

    dt
}
simulate_data.plcp_multi <- function(paras, n = 1) {
    simulate_data.plcp(paras[n,])
}

#' Perform a simulation study using a \code{study_parameters}-object
#'
#' @param object An object created by \code{\link{study_parameters}}.
#' @param nsim The number of simulations to run.
#' @param seed Currently ignored.
#' @param formula \code{lme4::lmer} formula(s) used to analyze the data, see \emph{Details}.
#' Should either be a character vector if one model is simulated.
#' It also possible to compare two models, e.g. a correct and a misspecified model,
#' by passing the two formulas as a named list, with the names "wrong" and "correct".
#' See \emph{Examples}.
#' @param satterthwaite Logical; if \code{TRUE} Satterthwaite's degrees of freedom
#' approximation will be used when computing \emph{p}-values. This is implemented using
#' the \code{lmerTest}-package. See \emph{Details}.
#' @param CI Logical; if \code{TRUE} coverage rates for confidence intervals
#' will be calculated. See \emph{Details}.
#' @param cores Number of CPU cores to use. Parallelization is done using
#' parallel::mclapply() which does not support \strong{Windows}. Use \code{cores = 1}
#' if you're on Windows.
#' @param progress \code{logical}; will display progress if \code{TRUE}. Currently
#' ignored on \emph{Windows}. Package \code{pbmclapply} is used to display progress.
#' \strong{N.B} using a progress bar will noticeably increase the simulation time.
#' @param batch_progress \code{logical}; if \code{TRUE} progress will be shown for
#' simulations with multiple setups.
#' @param save Logical; if \code{TRUE} results will be saved to your working directory
#' under the folder "save".
#' @param ... Optional arguments.
#'
#' @importFrom stats simulate as.formula confint pnorm pt qnorm qt rmultinom sd time vcov
#'
#' @details
#'
#' \strong{Model formula}
#'
#' The available model terms are:
#' \itemize{
#'  \item \code{y} the outcome vector
#'  \item \code{time} the time vector
#'  \item \code{treatment} treatment indicator (0 = "control", 1 = "treatment")
#'  \item \code{subject} subject-level id variable, from 1 to total number of subjects.
#'  \item \code{cluster} for three-level models; the cluster-level id variable,
#'  from 1 to the total number of clusters.
#' }
#'
#' See \emph{Examples} and the simulation-vignette for formula examples.
#'
#' \strong{Satterthwaite's approximation, and CI coverage}
#'
#' To decrease the simulation time Satterthwaite's dfs and the CIs' coverage rates
#' will only be calculated for the test of 'time:treatment'-interaction.
#'
#' Confidence intervals are both calculated using profile likelihood and by
#' the Wald approximation.
#'
#' @seealso \code{\link{summary.plcp_sim}}, \code{\link{simulate_data}}
#'
#' @examples
#' \dontrun{
#' # Two-level ---------------------------------------------------------------
#' p <- study_parameters(n1 = 11,
#'                       n2 = 25,
#'                       sigma_subject_intercept = 1.44,
#'                       sigma_subject_slope = 0.2,
#'                       sigma_error = 1.44,
#'                       cohend = 0.5)
#'
#' f <- "y ~ treatment * time + (1 + time | subject)"
#'
#'
#' res <- simulate(object = p,
#'                 nsim = 1000,
#'                 formula = f,
#'                 satterthwaite = TRUE,
#'                 progress = FALSE,
#'                 cores = 1,
#'                 save = FALSE)
#'
#' summary(res)
#'
#'
#' # Three-level (nested) ---------------------------------------------------------
#' p <- study_parameters(n1 = 10,
#'                       n2 = 20,
#'                       n3 = 4,
#'                       sigma_subject_intercept = 1.44,
#'                       icc_pre_cluster = 0,
#'                       sigma_subject_slope = 0.2,
#'                       icc_slope = 0.05,
#'                       sigma_error = 1.44,
#'                       cohend = 0.5)
#'
#' ## compare correct and miss-specified model
#' f <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (time | cluster)",
#'           "wrong" = "y ~ treatment * time + (1 + time | subject)")
#'
#' res <- simulate(object = p,
#'                 nsim = 1000,
#'                 formula = f,
#'                 satterthwaite = TRUE,
#'                 progress = FALSE,
#'                 cores = 1,
#'                 save = FALSE)
#'
#' summary(res)
#'
#'
#' # Partially nested design ------------------------------------------------
#' p <- study_parameters(n1 = 11,
#'                       n2 = 10,
#'                       n3 = 4,
#'                       sigma_subject_intercept = 1.44,
#'                       icc_pre_cluster = 0,
#'                       sigma_subject_slope = 0.2,
#'                       cor_subject = -0.5,
#'                       icc_slope = 0.05,
#'                       sigma_error = 1.44,
#'                       partially_nested = TRUE,
#'                       cohend = -0.5)
#'
#' f <- "y ~ treatment * time + (1 + time | subject) + (0 + treatment:time | cluster)"
#'
#' res <- simulate(object = p,
#'                 nsim = 1000,
#'                 formula = f,
#'                 satterthwaite = TRUE,
#'                 progress = FALSE,
#'                 cores = 4,
#'                 save = FALSE)
#'
#' summary(res)
#'
#' # Run multiple designs  --------------------------------------------------------
#' p <- study_parameters(n1 = 10,
#'                       n2 = 20,
#'                       n3 = c(2, 4, 6),
#'                       sigma_subject_intercept = 1.44,
#'                       icc_pre_cluster = 0,
#'                       sigma_subject_slope = 0.2,
#'                       icc_slope = 0.05,
#'                       sigma_error = 1.44,
#'                       cohend = 0.5)
#'
#' f <- list("correct" = "y ~ treatment * time + (1 + time | subject) + (time | cluster)",
#'           "wrong" = "y ~ treatment * time + (1 + time | subject)")
#'
#' res <- simulate(object = p,
#'                 nsim = 1000,
#'                 formula = f,
#'                 satterthwaite = TRUE,
#'                 progress = FALSE,
#'                 cores = 1,
#'                 save = FALSE)
#'
#' # Summarize 'time:treatment' results for n3 = c(2, 4, 6) for 'correct' model
#' summary(res, para = "time:treatment", type = "fixed", model = "correct")
#'
#' # Summarize cluster-level random slope  for n3 = c(2, 4, 6) for 'correct' model
#' summary(res, para = "cluster_slope", type = "random", model = "correct")
#' }
#' @export
simulate.plcp <- function(object,
                          nsim,
                          seed = NULL,
                          formula,
                          satterthwaite = FALSE,
                          CI = FALSE,
                          cores = 1,
                          progress = FALSE,
                          batch_progress = TRUE,
                          save = FALSE,
                          ...) {
    formula <- check_formula(formula)

    if (is.null(nrow(object))) {
        simulate.plcp_list(
            object = object,
            nsim = nsim,
            formula = formula,
            satterthwaite = satterthwaite,
            CI = CI,
            cores = cores,
            progress = progress,
            save = save,
            ...
        )
    } else if (nrow(object) > 1) {
        simulate.plcp_data_frame(
            object = object,
            nsim = nsim,
            formula = formula,
            satterthwaite = satterthwaite,
            CI = CI,
            cores = cores,
            progress = progress,
            batch_progress = batch_progress,
            save = save,
            ...
        )
    }
}


#' @rdname simulate.plcp
#' @export
simulate.plcp_multi <- function(object,
                                nsim,
                                seed = NULL,
                                formula,
                                satterthwaite = FALSE,
                                CI = FALSE,
                                cores = 1,
                                progress = FALSE,
                                batch_progress = TRUE,
                                save = FALSE,
                                ...) {
    simulate.plcp(
        object = object,
        nsim = nsim,
        seed = NULL,
        satterthwaite = satterthwaite,
        CI = CI,
        formula = formula,
        cores = cores,
        progress = progress,
        batch_progress = batch_progress,
        save = save,
        ...
    )
}
simulate.plcp_list <-
    function(object,
             nsim,
             seed = NULL,
             formula,
             satterthwaite = FALSE,
             CI = FALSE,
             cores = 1,
             progress = FALSE,
             save = FALSE,
             ...) {
        ptm <- proc.time()
        is_windows <- .Platform$OS.type == "windows"
        if (progress & !is_windows) {
            res <-
                pbmcapply::pbmclapply(
                    1:nsim,
                    simulate_,
                    paras = object,
                    formula = formula,
                    satterthwaite = satterthwaite,
                    CI = CI,
                    mc.cores = cores,
                    ...
                )
        } else {
            res <-
                parallel::mclapply(
                    1:nsim,
                    simulate_,
                    paras = object,
                    formula = formula,
                    satterthwaite = satterthwaite,
                    CI = CI,
                    mc.cores = cores,
                    ...
                )
        }
        time <- proc.time() - ptm


        out <-
            list(
                res = res,
                paras = object,
                nsim = nsim,
                time = time["elapsed"]
            )
        out <- munge_results(out)
        class(out) <- append(class(out), "plcp_sim")

        out
    }
simulate.plcp_data_frame <-
    function(object,
             nsim,
             seed = NULL,
             formula,
             satterthwaite,
             CI = FALSE,
             cores,
             progress,
             save = FALSE,
             batch_progress = TRUE) {
        if (save) {
            output_dir <- format(Sys.time(), "%Y%m%d_%H%M_%S")
            output_dir <- paste("save/", output_dir, sep = "")
            dir.create(output_dir)

            cat("Each batch is being saved to: ", output_dir, "/\n", sep = "")
        }

        res <- lapply(1:nrow(object), function(i) {
            if (batch_progress) {
                cat("\rBatch: ",
                    i,
                    "/",
                    nrow(object),
                    rep(" ", options()$width - 13))
                cat("\n")
            }

            x <- simulate.plcp_list(
                as.plcp_sim(object[i,]),
                nsim = nsim,
                formula = formula,
                satterthwaite = satterthwaite,
                CI = CI,
                cores = cores,
                progress = progress
            )

            if (save) {
                f <- paste(output_dir, "/sim", i, ".rds", sep = "")
                saveRDS(x, f)
            }
            x
        })

        class(res) <- append(class(res), "plcp_multi_sim")
        res
    }
simulate_ <- function(sim, paras, satterthwaite, CI, formula) {
    d <- simulate_data(paras)
    fit <- analyze_data(formula, d)
    res <- extract_results(fit, CI, satterthwaite = satterthwaite)
    res
}


# Checks ------------------------------------------------------------------

# checks
check_formula <- function(formula) {
    if (is.character(formula))
        formula <- list("correct" = formula)
    if (is.null(names(formula)))
        stop("Formula should be a named list (either 'correct' or 'wrong')")
    if (!all(names(formula) %in% c("correct", "wrong")))
        stop("Formula names must be either 'correct' or 'wrong'")

    if (length(names(formula)) > 1) {
        if (length(names(formula)) != length(unique(names(formula)))) {
            stop(
                "Both formulas can't have the same name. Allowed names are 'correct' or 'wrong'"
            )
        }
    }



    for (f in formula)
        check_formula_terms(f)

    formula
}
check_formula_terms <- function(f) {
    f <- as.formula(f)

    x <- all.vars(f)
    ind <- x %in% c("y", "treatment", "time", "subject", "cluster")
    wrong <- x[!ind]

    if (length(wrong) > 0) {
        stop(
            paste(
                wrong,
                "is not an allowed variable name. Allowed names are: 'y', 'treatment', 'time', 'subject', 'cluster'"
            )
        )
    }

    # test functions
    x <- all.vars(f, functions = TRUE)
    ind <- x %in% c("y", "treatment", "time", "subject", "cluster")
    x <- x[!ind]
    ind <- x %in% c("~", "+", "*", "(", "|", ":")
    x <- x[!ind]
    if (length(x) > 0) {
        stop(
            paste(
                "Functions are not yet supported in model formula. Problematic function was:",
                x
            )
        )
    }
}



# analyze -----------------------------------------------------------------

analyze_data <- function(formula, d) {

        fit <-
            lapply(formula, function(f)
                tryCatch(
                    #do.call(lme4::lmer, list(formula=f, data=d))
                    lme4::lmer(formula = as.formula(f), data = d)
                    )
                )



    fit
}

extract_results <- function(fit, CI = FALSE, satterthwaite = FALSE) {
    lapply(fit, extract_results_, CI = CI, satterthwaite = satterthwaite)
}
extract_random_effects <- function(fit) {
    x <- as.data.frame(lme4::VarCorr(fit))

    vcovs <- x[is.na(x$var2),]
    vcovs$parameter <- paste(vcovs$grp, vcovs$var1, sep = "_")
    vcovs$sdcor <- NULL

    correlations <-  x[!is.na(x$var2),]
    correlations$parameter <- paste(correlations$grp,
                                    correlations$var1,
                                    correlations$var2,
                                    sep = "_")
    correlations$vcov <- NULL

    names(correlations)[names(correlations) == "sdcor"] <- "vcov"

    rbind(vcovs, correlations)
}
add_p_value <- function(fit, satterthwaite) {
    tmp <- tryCatch(summary(fit))
    tmp <- tmp$coefficients

    ## satterthwaite
    if(satterthwaite) {

        ff <- rownames(tmp)
        ind <- which(ff %in%  c("treatment:time", "time:treatment"))

        L <- rep(0, length(ff))
        L[ind] <- 1

        satt <- tryCatch(lmerTest::calcSatterth(fit, L),  error = function(e) { NA })

        df <- rep(NA, length(ff))
        p <- rep(NA, length(ff))

        if(is.list(satt)) {
            df[ind] <- satt$denom
            p[ind] <- satt$pvalue
        }


        res <- list("df" = df,
                    "p" = p)
    } else {
        res <- list("df" = NA,
                    "p" = NA)
    }


    res

}
extract_results_ <- function(fit, CI, satterthwaite) {
    se <- sqrt(diag(vcov(fit)))
    tmp_p <- add_p_value(fit, satterthwaite)
    FE <- data.frame(
        "estimate" = lme4::fixef(fit),
        "se" = se,
        "pval" = tmp_p$p,
        "df" = tmp_p$df
    )


    rnames <- rownames(FE)
    parm_id <- which(c("time:treatment", "treatment:time") %in% rnames)
    CI_parm <- c("time:treatment", "treatment:time")[parm_id]
    rnames <- gsub("treatment:time", "time:treatment", rnames)
    rownames(FE) <- NULL
    FE <- cbind(data.frame(parameter = rnames), FE)

    if (CI) {
        CI <- tryCatch(stats::confint(fit, parm = CI_parm),
                       error = function(e) NA)
        CI_wald <-
            tryCatch(stats::confint(fit, method = "Wald", parm = CI_parm),
                     error = function(e) NA)
        FE[FE$parameter == "time:treatment", "CI_lwr"] <- CI[1]
        FE[FE$parameter == "time:treatment", "CI_upr"] <- CI[2]
        FE[FE$parameter == "time:treatment", "CI_wald_lwr"] <- CI_wald[1]
        FE[FE$parameter == "time:treatment", "CI_wald_upr"] <- CI_wald[2]

    }


    RE <- extract_random_effects(fit)

    conv <- is.null(fit@optinfo$conv$lme4$code)

    list("RE" = RE,
         "FE" = FE,
         "conv" = conv)
}





.rename_rr <- function(.x, match, new) {
    .x[.x$parameter %in% match, "parameter"] <- new

    .x
}
rename_random_effects <- function(.x) {
    .x <- .rename_rr(
        .x,
        match = c(
            "cluster_time",
            "cluster.1_time",
            "cluster_time:treatment",
            "cluster_treatment:time"
        ),
        new = "cluster_slope"
    )

    .x <- .rename_rr(.x,
                     match = "Residual_NA",
                     new = "error")

    .x <- .rename_rr(.x,
                     match = c("subject_(Intercept)"),
                     new = "subject_intercept")

    .x <- .rename_rr(.x,
                     match = c("cluster_(Intercept)"),
                     new = "cluster_intercept")

    .x <- .rename_rr(.x,
                     match = c("subject.1_time", "subject_time"),
                     new = "subject_slope")


    .x <- .rename_rr(.x,
                     match = "subject_(Intercept)_time",
                     new = "cor_subject")

    .x <- .rename_rr(.x,
                     match = c("cluster_(Intercept)_time",
                               "cluster_(Intercept)_treatment:time",
                               "cluster_(Intercept)_time:treatment"),
                     new = "cor_cluster")

}

munge_results <- function(res) {
    x <- res$res

    models <- names(x[[1]])
    RE <- lapply(models, munge_results_, x, "RE")
    FE <- lapply(models, munge_results_, x, "FE")

    convergence <- lapply(models, munge_results_, x, "conv")

    RE <- lapply(RE, rename_random_effects)

    x <- list()

    for (i in seq_along(models)) {
        x[[i]] <- list("RE" = RE[[i]],
                       "FE" = FE[[i]],
                       "convergence" = convergence[[i]])
    }

    names(x) <- models
    res$res <- x

    res
}
munge_results_ <- function(model, res, effect) {
    res <- lapply(seq_along(res), function(i) {
        x <- res[[i]][[model]][[effect]]

        x
    })
    res <- do.call(rbind, res)

    res
}



# Print -------------------------------------------------------------------
print_model <- function(i, x) {
    models <- names(x)
    cat("Model: ", models[i], "\n")
    cat("  Random effects\n\n")
    print.data.frame(x[[i]]$RE,
                     row.names = FALSE,
                     digits = 3,
                     quote = FALSE)
    cat("\n  Fixed effects\n\n")
    print.data.frame(x[[i]]$FE, row.names = FALSE, digits = 3)
    cat("---\n")

}
#' Print method for \code{summary.plcp_sim}-objects
#' @param x An object of class \code{plcp_sim_summary}
#' @param ... Optional arguments.
#' @method print plcp_sim_summary
#' @export
print.plcp_sim_summary <- function(x, ...) {
    res <- x
    #print
    x <- res$summary

    tot_n <- get_tot_n(res$paras)$total
    n3 <- get_n3(res$paras)
    lapply(seq_along(x), print_model, x)
    cat("Number of simulations:", res$nsim)
    cat("\nTime points (n1): ", res$paras$n1)
    cat(" | Subjects per cluster (n2): ", unlist(res$paras$n2))
    cat(" | Clusters per treatment (n3): ", unlist(n3[c("treatment", "control")]))
    cat("\nTotal number of subjects: ", tot_n)
    convergence <- lapply(x, function(d)
        d$convergence)
    convergence <- unlist(convergence)
    cat("\nConvergence: ",
        paste(convergence * 100, "%", collapse = ", "))
}

#' Print method for \code{simulate.plcp_multi}-objects
#' @param x An object created with \code{\link{simulate.plcp_multi}}
#' @param ... Optional arguments.
#' @importFrom utils object.size
#' @method print plcp_multi_sim
#' @export
print.plcp_multi_sim <- function(x, ...) {
    nsim <- x[[1]]$nsim
    nmulti <- length(x)
    cat("# ", nmulti, " x ", nsim, "grid of simulations")
    cat("\n")
    cat("# Object size:", format(object.size(x), units = "auto"))
}

#' Print method for \code{simulate.plcp}-objects
#' @param x An object created with \code{\link{simulate.plcp}}
#' @param ... Optional arguments.
#' @method print plcp_sim
#' @export
print.plcp_sim <- function(x, ...) {
    cat(
        "# A 'plcp_sim'-object containing",
        x$nsim,
        "simulations.",
        "Use summary() to view results."
    )
    cat("\n")
    cat("# Object size:", format(object.size(x), units = "auto"))
    invisible(x)
}


# #' @export
# p_sim <- function(res) {
#     res <- lapply(res$res, summary_.plcp_sim, paras = res$paras)
#
#
# }

#' Summarize the results from a simulation of a single study design-object
#' @param object A \code{simulate.plcp}-object
#'
#' @param ... Currently not used
#'
#' @return Object with class \code{plcp_sim_summary}. It contains
#' the following output:
#' \itemize{
#'  \item \code{parameter} is the name of the coefficient
#'  \item \code{M_est} is the mean of the estimates taken over all the simulations.
#'  \item \code{M_se} is the mean estimated standard error taken over all the simulations.
#'  \item \code{SD_est} is the empirical standard error; i.e. the standard
#'  deviation of the distribution of the generated estimates.
#'  \item \code{power} is the empirical power of the Wald Z test, i.e. the proportion
#'  of simulated p-values < alpha.
#'  \item \code{power_satt} is the empirical power of the Wald \emph{t} test using
#'   Satterthwaite's degree of freedom approximation.
#'  \item \code{satt_NA} is the proportion of Satterthwaite's approximations that failed.
#'  \item \code{prop_zero} is the proportion of the simulated estimates that
#'  are zero; only shown for random effects.
#' }
#'
#'
#' @method summary plcp_sim
#' @export
summary.plcp_sim <- function(object, ...) {
    res <- object
    x <- lapply(res$res, summary_.plcp_sim, paras = res$paras)
    x <- list(summary = x,
              nsim = res$nsim,
              paras = res$paras)
    class(x) <- append("plcp_sim_summary", class(x))
    x
}
summary_.plcp_sim  <- function(res, paras) {
    RE_params <-
        data.frame(
            parameter = c(
                "subject_intercept",
                "subject_slope",
                "cluster_intercept",
                "cluster_slope",
                "error",
                "cor_subject",
                "cor_cluster"
            ),
            theta = c(
                paras$sigma_subject_intercept ^ 2,
                paras$sigma_subject_slope ^ 2,
                paras$sigma_cluster_intercept ^ 2,
                paras$sigma_cluster_slope ^ 2,
                paras$sigma_error ^ 2,
                paras$cor_subject,
                paras$cor_cluster
            )
        )


    RE <- lapply(unique(res$RE$parameter), function(i, .d, theta) {
        vcov <-  .d[.d$parameter == i, "vcov"]

        para <- i
        theta <- theta[theta$parameter == i, "theta"]
        res <- data.frame(
            parameter = para,
            M_est = mean(vcov),
            theta = theta,
            prop_zero = mean(abs(vcov - 0) < .Machine$double.eps ^
                                 0.5)
        )
        if (para %in% c("cor_subject", "cor_cluster")) {
            res$prop_zero <- mean(abs(abs(vcov) - 1) < .Machine$double.eps ^ 0.5)
        }
        res

    }, .d = res$RE, theta = RE_params)

    RE <- do.call(rbind, RE)

    # df
    false_conv <-
        res$RE[res$RE$parameter == "cluster_slope", "vcov"]
    false_conv <-
        mean(abs(false_conv - 0 < .Machine$double.eps ^ 0.5))


    df <- ifelse(false_conv, get_tot_n(paras) - 2, get_n3(paras) - 2)


    theta = list(
        "(Intercept)" = paras$fixed_intercept,
        "treatment" = 0,
        "time" = paras$fixed_slope,
        "time:treatment" = get_slope_diff(paras) / paras$T_end
    )
    FE <- lapply(unique(res$FE$parameter), function(i, .d, theta) {
        se <- .d[.d$parameter == i, "se"]
        estimate <- .d[.d$parameter == i, "estimate"]
        pval <- .d[.d$parameter == i, "pval"]
        df <- .d[.d$parameter == i, "df"]

        para <- i
        theta <- theta[[i]]
        data.frame(
            parameter = i,
            M_est = mean(estimate),
            M_se = mean(se),
            SD_est = sd(estimate),
            Power = mean(get_cover(estimate, se)),
            Power_satt = mean(pval < 0.05, na.rm = TRUE),
            Satt_NA = mean(is.na(pval))
        )
    }, .d = res$FE, theta = theta)

    FE <- do.call(rbind, FE)

    ## Fix this, remove dplyr
    if ("CI_lwr" %in% colnames(res$FE)) {
        diff <- get_slope_diff(paras) / paras$T_end

        CI_cov <- lapply(unique(res$FE$parameter), function(i, .d, diff) {
            CI_lwr <- .d[.d$parameter == i, "CI_lwr"]
            CI_upr <- .d[.d$parameter == i, "CI_upr"]
            CI_wald_lwr <- .d[.d$parameter == i, "CI_wald_lwr"]
            CI_wald_upr <- .d[.d$parameter == i, "CI_wald_upr"]

            para <- i
            data.frame(
                parameter = i,
                CI_cover =  mean(CI_lwr < diff &
                                     CI_upr > diff, na.rm=TRUE),
                CI_Wald_cover = mean(CI_wald_lwr < diff &
                                         CI_wald_upr > diff, na.rm=TRUE))

        }, .d = res$FE, diff = diff)

        CI_cov <- do.call(rbind, CI_cov)

        FE$CI_Cover <- CI_cov$CI_cover
        FE$CI_Wald_cover <- CI_cov$CI_Wald_cover
    }


    FE_eff <- c("(Intercept)",
                "treatment",
                "time",
                "time:treatment")

    FE_eff <- FE_eff[FE_eff %in% FE$parameter]
    ind <- vapply(FE_eff, function(i)
                        which(FE$parameter == i), numeric(1))
    FE <- FE[ind,]

    convergence <- mean(res$convergence[[2]])

    list("RE" = RE,
         "FE" = FE,
         "convergence" = convergence)
}


# multi-sim ---------------------------------------------------------------


as.plcp_sim <- function(object) {
    UseMethod("as.plcp_sim")
}

as.plcp_sim.data.frame <- function(object) {
    object <- as.list(object)
    class(object) <- append("plcp", class(object))

    object
}


#' Summarize simulations based on a combination of multiple parameter values
#'
#' @param object A multiple simulation object created with
#' \code{\link{simulate.plcp_multi}}
#' @param para The name of the fixed or random effect that should be summarized.
#' @param type Specifies what type of effect \code{para} is; can be either
#' "random" or "fixed".
#' @param model Specifies which model should be summarized. Used when models are
#' simulated from formulas named "correct" and "wrong".
#' @param ... Optional arguments.
#'
#' @method summary plcp_multi_sim
#'
#' @return A \code{data.frame} with class \code{plcp_multi_sim_summary}. Each row
#' summarizes one of the parameter combinations used in the simulation. It contains
#' the following columns:
#' \itemize{
#'  \item \code{parameter} is the name of the coefficient
#'  \item \code{M_est} is the mean of the estimates taken over all the simulations.
#'  \item \code{theta} is the population parameter values specified with \code{study_parameters}
#'  \item \code{est_rel_bias} is the relative bias of the estimate
#'  \item \code{M_se} is the mean estimated standard error taken over all the simulations.
#'  \item \code{SD_est} is the empirical standard error; i.e. the standard
#'  deviation of the distribution of the generated estimates
#'  \item \code{se_rel_bias} is the relative bias of the standard error.
#'  \item \code{power} is the empirical power of the Wald Z test, i.e. the proportion
#'  of simulated p-values < alpha
#'  \item \code{power_satt} is the empirical power of the Wald \emph{t} test using
#'   Satterthwaite's degree of freedom approximation
#'  \item \code{satt_NA} is the proportion of Satterthwaite's approximations that failed.
#'  \item \code{prop_zero} is the proportion of the simulated estimates that
#'  are zero; only shown for random effects
#' }
#'
#' @export
#'
summary.plcp_multi_sim <- function(object,
                                   para = "time:treatment",
                                   type = "fixed",
                                   model = "correct", ...) {
    res <- object
    if (!model %in% (names(res[[1]][[1]])))
        stop("No model named: ", model)
    if (!type %in% c("fixed", "random"))
        stop("'type' should be either 'fixed' or 'random'")
    if (type == "fixed" &
        !para %in% c("intercept", "time", "treatment", "time:treatment")) {
        stop("Para should be one of: 'intercept' 'time', 'treatment', 'time:treatment'")
    }
    if (type == "fixed") {
        out <- lapply(res, summary_fixed.plcp_multi_sim, para, model)
    } else {
        out <- lapply(res, summary_random.plcp_multi_sim, para, model)
    }

    out <- do.call(rbind, out)

    class(out) <- append("plcp_multi_sim_summary", class(out))
    attr(out, "type") <- type
    attr(out, "model") <- model
    attr(out, "nsim") <- res[[1]]$nsim
    out
}

#' @method print plcp_multi_sim_summary
#' @export
print.plcp_multi_sim_summary <- function(x, ...) {
    model <- attr(x, "model")
    type <- attr(x, "type")
    nsim <- attr(x, "nsim")
    x <- as.data.frame(x)
    cat("Model: ", model, "| Type:", type, "\n")
    cat("---\n")
    print(x, digits = 3)
    cat("---\n")
    cat("nsim: ", nsim)
    invisible(x)
}

## random effect
summary_fixed.plcp_multi_sim <- function(res, para, model) {

    theta <- switch(
        para,
        "intercept" = res$paras$fixed_intercept,
        "treatment" = 0,
        "time" = res$paras$fixed_slope,
        "time:treatment" = get_slope_diff(res$paras) / res$paras$T_end
    )

    out <- res$res[[model]]$FE
    out <- out[out$parameter == para,]


    est <- mean(out$estimate)
    se_est <- mean(out$se)
    se_hat <- sd(out$estimate)

    power <- mean(get_cover(out$estimate, out$se))
    out <- with(
        out,
        data.frame(
            parameter = para,
            M_est = est,
            theta = theta,
            est_rel_bias = ifelse(theta == 0,
                                  est - theta,
                                  (est - theta) / theta),
            M_se = se_est,
            SD_est = se_hat,
            se_rel_bias = (se_est - se_hat) / se_hat,
            Power = power,
            Power_satt = mean(out$pval < 0.05, na.rm = TRUE),
            Satt_NA = mean(is.na(out$pval))

        )
    )

    out
}
summary_random.plcp_multi_sim <- function(res, para, model) {
    tmp <- res$res[[model]]$RE
    if (!para %in% unique(tmp$parameter)) {
        stop(
            "No random effect named: '",
            para,
            "'.\n Available random effects are: \n  ",
            paste(unique(tmp$parameter), collapse = ",\n  ")
        )
    }

    theta <- switch(
        para,
        "subject_intercept" = res$paras$sigma_subject_intercept,
        "subject_slope" = res$paras$sigma_subject_slope,
        "cluster_slope" = res$paras$sigma_cluster_slope,
        "error" = res$paras$sigma_error
    )

    out <- res$res[[model]]$RE
    out <- out[out$parameter == para,]

    est <- mean(out$vcov)
    out <- data.frame(
        parameter = para,
        M_est = est,
        theta = theta ^ 2,
        est_rel_bias = (est - theta ^ 2) / theta ^ 2,
        prop_zero = mean(abs(out$vcov - 0 < .Machine$double.eps ^
                                 0.5))
    )


    #out <- cbind(get_rel_bias(res$para), out)
    #out <- dplyr::select(out, -approx_rel_bias)

    out

}

# power -------------------------------------------------------------------
get_cover <- function(est, se) {
    # 95 % CI cover
    abs(est) - qnorm(0.975) * se > 0
}

get_p_val_df <- function(t, df, parameter) {
    if (unique(parameter) %in% c("time:treatment", "time")) {
        p <- (1 - pt(abs(t), df)) * 2
    } else
        p <- NA

    p
}
