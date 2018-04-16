#' Generate a data set using a \code{study_parameters}-object
#'
#' @param paras An object created by \code{\link{study_parameters}}
#' @param n Optional; specifies which row \code{n} should be used if \code{object}
#'  is a \code{data.frame} containing multiple setups.
#'
#' @return A \code{data.frame} with the simulated data in long form. With the following columns:
#' \itemize{
#'  \item \code{y} the outcome vector, with missing values as NA
#'  \item \code{y_c} the outcome vector, without missing values removed.
#'  \item \code{time} the time vector
#'  \item \code{treatment} treatment indicator (0 = "control", 1 = "treatment")
#'  \item \code{subject} subject-level id variable, from 1 to total number of subjects.
#'  \item \code{cluster} for three-level models; the cluster-level id variable,
#'  from 1 to the total number of clusters.
#' }
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

create_dropout_indicator <- function(paras) {
    ind <- dropout_process(unlist(paras$dropout), paras)
    ind <- ind$missing

    ind
}
add_NA_values_from_indicator <- function(d, missing) {
    d$miss <- missing
    d$y <- ifelse(d$miss == 1, NA, d$y)
    d <- d[, c("y",
               "y_c",
                grep("^y$|^y_c$", colnames(d), value = TRUE, invert = TRUE))]
    d
}
#' @rdname simulate_data
#' @export
simulate_data.plcp <- function(paras, n = NULL) {
    if (is.data.frame(paras))
        paras <- as.list(paras)
    if(is.null(paras$prepared)) {
        tmp <- prepare_paras(paras)
    } else tmp <- paras
    paras <- tmp$control
    paras_tx <- tmp$treatment

    slope_diff <- get_slope_diff(paras) / paras$T_end
    paras$effect_size <- NULL
    paras_tx$effect_size <- NULL

    paras_tx$fixed_slope <- paras_tx$fixed_slope + slope_diff

    paras_tx$partially_nested <- NULL
    paras$partially_nested <- NULL
    paras$allocation_ratio <- NULL
    paras_tx$allocation_ratio <- NULL

    # replace NA
    paras[is.na(paras)] <- 0
    paras_tx[is.na(paras_tx)] <- 0

    d_tx <- simulate_3lvl_data(paras_tx)
    d_c <- simulate_3lvl_data(paras)

    # drop outs
    if (is.list(paras$dropout) |
        is.function(paras$dropout) |
        length(paras$dropout) > 1) {

        miss_c <- create_dropout_indicator(paras)
        d_c <- add_NA_values_from_indicator(d_c, miss_c)
    }
    if (is.list(paras_tx$dropout) |
        is.function(paras_tx$dropout) |
        length(paras_tx$dropout) > 1) {
        miss_tx <- create_dropout_indicator(paras_tx)
        d_tx <- add_NA_values_from_indicator(d_tx, miss_tx)
    }

    # combine
    d_tx$treatment <- 1
    d_c$treatment <- 0
    d_c$subject <- d_c$subject + max(d_tx$subject)
    d_tx$cluster <- d_tx$cluster + max(d_c$cluster)
    dt <- rbind(d_c, d_tx)

    dt
}
#' @rdname simulate_data
#' @export
simulate_data.plcp_multi <- function(paras, n = 1) {
    simulate_data.plcp(as.plcp(paras[n,]))
}

#' Perform a simulation study using a \code{study_parameters}-object
#'
#' @param object An object created by \code{\link{study_parameters}}.
#' @param nsim The number of simulations to run.
#' @param seed Currently ignored.
#' @param formula \code{lme4::lmer} formula(s) used to analyze the data, see \emph{Details}.
#' Should either be a character vector or named list if one model is simulated.
#' It also possible to compare two models, e.g. a correct and a misspecified model,
#' by passing the two formulas as a named list, with the names "wrong" and "correct".
#' See \emph{Examples}. If \code{NULL} the formula is made automatically,
#' using \code{\link{create_lmer_formula}}, which does not support objects with multiple
#' simulation setups.
#' @param satterthwaite Logical; if \code{TRUE} Satterthwaite's degrees of freedom
#' approximation will be used when computing \emph{p}-values. This is implemented using
#' the \code{lmerTest}-package. See \emph{Details}.
#' @param CI Logical; if \code{TRUE} coverage rates for confidence intervals
#' will be calculated. See \emph{Details}.
#' @param cores Number of CPU cores to use. If called from a GUI environment (e.g. RStudio) or
#' a computer running Microsoft Windows, PSOCK clusters will be used. If called from a
#' non-interactive Unix environment forking is utilized.
#' @param progress \code{logical}; will display progress if \code{TRUE}. Currently
#' ignored on \emph{Windows}. Package \code{pbmclapply} is used to display progress,
#' which relies on forking. \strong{N.B} using a progress bar will noticeably
#' increase the simulation time, due to the added overhead.
#' @param batch_progress \code{logical}; if \code{TRUE} progress will be shown for
#' simulations with multiple setups.
#' @param ... Optional arguments, see \emph{Saving} in \emph{Details} section.
#'
#' @importFrom stats simulate as.formula confint pnorm pt qnorm qt rmultinom sd time vcov reshape
#'
#' @details
#'
#' \strong{Model formula}
#'
#' The available model terms are:
#' \itemize{
#'  \item \code{y} the outcome vector, with potential missing data.
#'  \item \code{y_c} the complete version of \code{y}, before dropout was simulated.
#'  \item \code{time} the time vector.
#'  \item \code{treatment} treatment indicator (0 = "control", 1 = "treatment").
#'  \item \code{subject} subject-level id variable, from 1 to total number of subjects.
#'  \item \code{cluster} for three-level models; the cluster-level id variable,
#'  from 1 to the total number of clusters.
#' }
#'
#' Currently, the models are assumed to be of the form "\code{y ~ time * treatment + (time | subject) + (time | clusters)}",
#' i.e., "\code{subject}" and "\code{cluster}" are included as random effects.
#' The \code{\link{summary.plcp_sim}} method will not
#' work if you include parameters that are not defined by the study_parameters-object,
#' e.g. "\code{y ~ time * cluster + (1 | subject)}".
#'
#' See \emph{Examples} and the simulation-vignette for formula examples. For
#' \code{object}s that contain a single study setup, then the lmer formula
#' can be created automatically using \code{\link{create_lmer_formula}}.
#'
#' \strong{Satterthwaite's approximation, and CI coverage}
#'
#' To decrease the simulation time Satterthwaite's \emph{dfs} and the CIs' coverage rates
#' will only be calculated for the test of 'time:treatment'-interaction.
#'
#' Confidence intervals are both calculated using profile likelihood and by
#' the Wald approximation.
#'
#' \strong{Saving intermediate results for multi-sims}
#'
#' Objects with multi-sims can be save after each batch is finished. This is highly
#' recommended when many designs are simulated. The following additional arguments
#' control saving behavior:
#'
#' \itemize{
#'  \item \code{'save'}, \code{logical}, if \code{TRUE} each batch is saved as a
#'  \code{RDS}-file. Results are saved in your working directory, in the directory
#'  specified by \code{save_folder}.
#'  \item \code{'save_folder'} a \code{character} indicating the folder name. Default is \code{'save'}.
#'  \item \code{'save_folder_create'},  \code{logical}, if \code{TRUE} then \code{save_folder}
#'  will be created if it does not exist in your working directory.
#' }
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
#'                       effect_size = cohend(0.5))
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
#'                       effect_size = cohend(0.5))
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
#'                       effect_size = cohend(-0.5))
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
#'                       effect_size = cohend(0.5))
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
                          formula = NULL,
                          satterthwaite = FALSE,
                          CI = FALSE,
                          cores = 1,
                          progress = FALSE,
                          batch_progress = TRUE,
                          ...) {

    if(is.null(formula)) formula <- create_lmer_formula(object)
    formula <- check_formula(formula)
    if(satterthwaite) check_installed("lmerTest")

    if (is.null(nrow(object))) {
        simulate.plcp_list(
            object = object,
            nsim = nsim,
            formula = formula,
            satterthwaite = satterthwaite,
            CI = CI,
            cores = cores,
            progress = progress,
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
            ...
        )
    }
}


#' @rdname simulate.plcp
#' @export
simulate.plcp_multi <- function(object,
                                nsim,
                                seed = NULL,
                                formula = NULL,
                                satterthwaite = FALSE,
                                CI = FALSE,
                                cores = 1,
                                progress = FALSE,
                                batch_progress = TRUE,
                                ...) {

    if(is.null(formula)) stop("Argument 'formula' is missing. Automatic formula creation is not yet supported for multi-sim objects")
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
             progress = TRUE,
             save = FALSE,
             save_folder = NULL,
             save_folder_create = NULL,
             cl = NULL,
             ...) {
        ptm <- proc.time()
        is_windows <- .Platform$OS.type == "windows"

        if (progress & !is_windows) {
            check_installed("pbmcapply")
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
        } else if (cores > 1) {
            if(.Platform$OS.type == "unix" && !interactive()) {
                res <- parallel::mclapply(X = 1:nsim,
                                 FUN = simulate_,
                                 paras = object,
                                 formula = formula,
                                 satterthwaite = satterthwaite,
                                 CI = CI,
                                 mc.cores = cores,
                                 ...)
            } else {
                if(is.null(cl)) {
                    cl <- parallel::makeCluster(min(cores, nsim))
                    on.exit(parallel::stopCluster(cl))
                }
                parallel::clusterEvalQ(cl, expr =
                                           suppressPackageStartupMessages(require(powerlmm, quietly = TRUE)))
                clust_call <- (function(object, satterthwaite, formula, CI) {
                    function(i) {
                        simulate_(i, paras = object,
                                  satterthwaite = satterthwaite,
                                  formula = formula,
                                  CI = CI)
                    }
                })(object, satterthwaite, formula, CI)

                res <- parLapply(cl,
                                 X = 1:nsim,
                                 clust_call
                                )
            }
        } else {
            res <- lapply(1:nsim,
                             FUN = simulate_,
                             paras = object,
                             formula = formula,
                             satterthwaite = satterthwaite,
                             CI = CI,
                             ...)

        }
        time <- proc.time() - ptm


        out <-
            list(
                res = res,
                paras = object,
                nsim = nsim,
                time = time["elapsed"],
                formula = formula
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
             batch_progress,
             save = FALSE,
             save_folder = "save",
             save_folder_create = FALSE) {

        if (save) {
            if(!dir.exists(save_folder)) {
                if(save_folder_create) {
                    message("Creating save_folder: ", file.path(save_folder))
                    dir.create(file.path(save_folder))
                } else {
                    stop("Directory '", save_folder, "' does not exist.")
                }
            }

            output_dir <- format(Sys.time(), "%Y%m%d_%H%M_%S")
            output_dir <- file.path(save_folder, output_dir)
            dir.create(output_dir)
            if(!dir.exists(output_dir)) stop("Could not create save directory.")

            cat("Each batch is being saved to: ", output_dir, "\n", sep = "")
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
            if(cores > 1 && interactive()) {
                cl <- parallel::makeCluster(cores)
                on.exit(parallel::stopCluster(cl))
            } else cl <- NULL
            x <- simulate.plcp_list(
                as.plcp(object[i,]),
                nsim = nsim,
                formula = formula,
                satterthwaite = satterthwaite,
                CI = CI,
                cores = cores,
                progress = progress,
                cl = cl
            )

            if (save) {
                fname <- paste("sim", i, ".rds", sep ="")
                f <- file.path(output_dir, fname)
                saveRDS(x, f)
            }
            x
        })
        if(save) saveRDS(object, file.path(output_dir, "paras.rds"))
        class(res) <- append(class(res), "plcp_multi_sim")
        res
    }
simulate_ <- function(sim, paras, satterthwaite, CI, formula) {
    prepped <- prepare_paras(paras)

    d <- simulate_data(prepped)
    #saveRDS(d, file = paste0("/tmp/R/sim",sim, ".rds"))
    tot_n <- length(unique(d[d$time == 0,]$subject))
    fit <- analyze_data(formula, d)
    res <- extract_results(fit = fit,
                           CI = CI,
                           satterthwaite = satterthwaite,
                           df_bw = get_balanced_df(prepped),
                           tot_n = tot_n)

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
    ind <- x %in% c("y", "y_c", "treatment", "time", "subject", "cluster")
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
    ind <- x %in% c("y", "y_c", "treatment", "time", "subject", "cluster")
    x <- x[!ind]
    ind <- x %in% c("~", "+", "*", "(", "|", "||", ":")
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

extract_results <- function(fit, CI = FALSE, satterthwaite = FALSE, df_bw, tot_n) {
    lapply(fit, extract_results_,
           CI = CI,
           satterthwaite = satterthwaite,
           df_bw = df_bw,
           tot_n = tot_n)
}
extract_random_effects <- function(fit) {
    x <- as.data.frame(lme4::VarCorr(fit))
    vcovs <- x[is.na(x$var2),]
    vcovs$parameter <- paste(vcovs$grp, vcovs$var1, sep = "_")
    vcovs$sdcor <- NULL

    correlations <-  x[!is.na(x$var2), ]
    correlations$parameter <- paste(correlations$grp,
                                    correlations$var1,
                                    correlations$var2,
                                    sep = "_")
    correlations$vcov <- NULL

    names(correlations)[names(correlations) == "sdcor"] <- "vcov"

    rbind(vcovs, correlations)
}

#' @importFrom utils packageVersion
add_p_value <- function(fit, satterthwaite, df_bw = NULL) {
    tmp <- tryCatch(summary(fit))
    tmp <- tmp$coefficients

    ff <- rownames(tmp)
    # satterth dfs only for time:treatment
    satterth_term <- any(ff %in%  c("treatment:time", "time:treatment"))
    ## satterthwaite
    if(satterthwaite & satterth_term) {
        ind <- which(ff %in%  c("treatment:time", "time:treatment"))

        L <- rep(0, length(ff))
        L[ind] <- 1

        # calcSatterth will be deprecated
        if(packageVersion("lmerTest") >= "3.0.0") {
            satt <- suppressMessages(tryCatch(lmerTest::contest(fit, L),
                                              error = function(e) { NA }))
        } else {
            satt <- suppressMessages(tryCatch(lmerTest::calcSatterth(fit, L),
                                              error = function(e) { NA }))
        }

        df <- rep(NA, length(ff))
        p <- rep(NA, length(ff))

        if(inherits(satt, "list")) {
            df[ind] <- satt$denom
            p[ind] <- satt$pvalue
        } else if(inherits(satt, "data.frame")) {
            df[ind] <- satt$DenDF
            p[ind] <- satt$`Pr(>F)`
        }
        if(is.na(p[ind])) {
            tval <- tmp[ind, "t value"]
            p[ind] <- 2*(1 - pt(abs(tval), df = df_bw))
        }

        res <- list("df" = df,
                    "p" = p)
    } else {
        res <- list("df" = NA,
                    "p" = NA)
    }


    res

}
fix_sath_NA_pval <- function(x, df) {
    ind <- is.na(x$pval)
    tmp <- x[ind, ]
    t <- with(tmp, estimate/se)
    pval <- 2*(1 - pt(abs(t), df = df))

    x[ind, "pval"] <- pval

    x
}
extract_results_ <- function(fit, CI, satterthwaite, df_bw, tot_n) {
    se <- sqrt(diag(vcov(fit)))
    tmp_p <- add_p_value(fit = fit,
                         satterthwaite = satterthwaite,
                         df_bw = df_bw)
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
    FE <- cbind(data.frame(parameter = rnames,
                           stringsAsFactors = FALSE),
                FE)

    FE[FE$parameter == "time:treatment", "df_bw"] <- df_bw

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
         "tot_n" = tot_n,
         "conv" = conv)
}

rename_rr_ <- function(.x, match, new) {
    .x[.x$parameter %in% match, "parameter"] <- new

    .x
}
rename_random_effects <- function(.x) {
    .x <- rename_rr_(
        .x,
        match = c(
            "cluster_time",
            "cluster.1_time",
            "cluster_time:treatment",
            "cluster_treatment:time",
            "cluster.1_treatment:time"
        ),
        new = "cluster_slope"
    )

    .x <- rename_rr_(.x,
                     match = "Residual_NA",
                     new = "error")

    .x <- rename_rr_(.x,
                     match = c("subject_(Intercept)"),
                     new = "subject_intercept")

    .x <- rename_rr_(.x,
                     match = c("cluster_(Intercept)",
                               "cluster_treatment"),
                     new = "cluster_intercept")

    .x <- rename_rr_(.x,
                     match = c("subject.1_time", "subject_time"),
                     new = "subject_slope")


    .x <- rename_rr_(.x,
                     match = "subject_(Intercept)_time",
                     new = "cor_subject")

    .x <- rename_rr_(.x,
                     match = c("cluster_(Intercept)_time",
                               "cluster_(Intercept)_treatment:time",
                               "cluster_(Intercept)_time:treatment",
                               "cluster_treatment_treatment:time"),
                     new = "cor_cluster")

}
order_model_lists <- function(models, RE, FE, tot_n, convergence) {
    x <- list()
    for (i in seq_along(models)) {
        x[[i]] <- list("RE" = RE[[i]],
                       "FE" = FE[[i]],
                       "tot_n" = tot_n[[i]],
                       "convergence" = convergence[[i]])
    }
    names(x) <- models

    x
}
munge_results <- function(res) {
    x <- res$res

    models <- names(x[[1]])
    RE <- lapply(models, munge_results_, x, "RE")
    FE <- lapply(models, munge_results_, x, "FE")
    convergence <- lapply(models, munge_results_, x, "conv")
    tot_n <- lapply(models, munge_results_, x, "tot_n")
    RE <- lapply(RE, rename_random_effects)

    x <- order_model_lists(models,
                           RE = RE,
                           FE = FE,
                           tot_n = tot_n,
                           convergence = convergence)

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
print_model <- function(i, x, digits) {
    models <- names(x)
    cat("Model: ", models[i], "\n")
    cat("  Random effects\n\n")
    print.data.frame(x[[i]]$RE,
                     row.names = FALSE,
                     digits = digits,
                     quote = FALSE)
    cat("\n  Fixed effects\n\n")
    print.data.frame(x[[i]]$FE, row.names = FALSE, digits = digits)
    cat("---\n")

}

#' @importFrom stats median
tot_n_ <- function(x) {

    x <- as.numeric(x)
    if(length(unique(x)) > 1) {
        x_mean <- mean(x)
        x_median <- median(x)
        x_sd <- sd(x)

        x <- paste(round(x_mean, 1),
                   "(mean)",
                   round(x_median, 1),
                   "(median)",
                   round(x_sd, 2),
                   "(SD)")
    }

    x
}

print_test_NA <- function(i, x) {
    mod_name <- names(x[i])
    Satt_NA <- x[[i]]$Satt_NA
    CI_NA <- x[[i]]$CI_NA
    non_convergence <- 1 - x[[i]]$convergence

    mess <- NULL
    if(Satt_NA > 0) {
        message("[Model: ", mod_name, "] ", round(Satt_NA, 4) * 100, "% of the Satterthwaite calculations failed")
    }
    if(CI_NA > 0) {
        message("[Model: ", mod_name, "] ", round(CI_NA, 4) * 100, "% of the profile likelihood CIs failed")
    }
    if(non_convergence > 0) {
        message("[Model: ", mod_name, "] ", round(non_convergence, 4) * 100, "% of the models threw convergence warnings")
    }
}

#' Print method for \code{summary.plcp_sim}-objects
#' @param x An object of class \code{plcp_sim_summary}
#' @param digits number of significant digits.
#' @param ... Optional arguments.
#' @method print plcp_sim_summary
#' @export
print.plcp_sim_summary <- function(x, digits = 2, ...) {
    res <- x
    #print
    x <- res$summary

    tot_n <- tot_n_(res$tot_n)
    if(length(unique(tot_n)) == 1) tot_n <- unique(tot_n)
    n3 <- get_n3(res$paras)
    n2 <- get_n2(res$paras)
    if(is.per_treatment(res$paras$n2)) {
        n2$treatment <- deparse_n2(n2$treatment)
        n2$control <- deparse_n2(n2$control)
        n2 <- unlist(n2)
    } else n2 <- deparse_n2(n2$treatment)


    if(is.unequal_clusters(res$paras$n2)) {
        n2_lab <- "\nSubjects per cluster: "
        n2 <- print_per_treatment(prepare_print_n2(res$paras),
                                  hanging = 23,
                                  n2 = TRUE)
    } else {
        n2_lab <- "\nSubjects per cluster (n2 x n3): "
        n2 <- print_per_treatment(prepare_print_n2(res$paras),
                                  hanging = 33,
                                  n2 = TRUE)
    }
    lapply(seq_along(x), print_model, x, digits = digits)
    cat("Number of simulations:", res$nsim, " | alpha: ", res$alpha)
    cat("\nTime points (n1): ", res$paras$n1)
    cat(n2_lab, n2)
    cat("\nTotal number of subjects: ", tot_n, "\n")
    lapply(seq_along(x), print_test_NA, x = x)
    if(any(x$correct$RE$is_NA > 0)) warning("Some estimated random effects were removed due to being NA.", call. = FALSE)
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
    cat("# lmer formula(s) used: \n")
    cat_formulas(x[[1]]$formula)
    cat("# Object size:", format(object.size(x), units = "auto"))
    cat("\n")
}

cat_formulas <- function(f) {
    nams <- names(f)
    names
    cat(paste("#    '", nams,"': ", f, sep = ""), sep = "\n")
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
    cat("# lmer formula(s) used: \n")
    cat_formulas(x$formula)
    cat("# Object size:", format(object.size(x), units = "auto"))
    cat("\n")
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
#' @param alpha Indicates the significance level. Default is 0.05 (two-tailed),
#' one-tailed tests are not yet implemented.
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
summary.plcp_sim <- function(object, alpha = 0.05, ...) {
    res <- object
    x <- lapply(res$res, summary_.plcp_sim,
                paras = res$paras,
                alpha = alpha)
    x <- list(summary = x,
              nsim = res$nsim,
              paras = res$paras,
              tot_n = x[[1]]$tot_n,
              alpha = alpha)
    class(x) <- append("plcp_sim_summary", class(x))
    x
}
summarize_RE <- function(res, theta) {
    d <- res$RE
    parms <- unique(d$parameter)
    tmp <- vector("list", length(parms))
    for(i in seq_along(parms)) {
        para <- parms[[i]]
        vcov <- d[d$parameter == para, "vcov"]
        theta_i <- theta[theta$parameter == para, "theta"]
        theta_i[is.na(theta_i)] <- 0 # for when para is NA
        if(length(theta_i) == 0) theta_i <- 0 # for when para does not exist
        est <- mean(vcov, na.rm=TRUE)
        res <- data.frame(
            parameter = para,
            M_est = est,
            theta = theta_i,
            est_rel_bias = get_RB(est, theta_i),
            prop_zero = mean(is_approx(vcov, 0), na.rm=TRUE),
            is_NA = mean(is.na(vcov))
        )
        if (para %in% c("cor_subject", "cor_cluster")) {
            res$prop_zero <- mean(is_approx(abs(vcov), 1), na.rm=TRUE)
        }
        tmp[[i]] <- res
    }
    RE <- do.call(rbind, tmp)

    RE
}
summarize_FE <- function(res, theta, alpha) {
    d <- res$FE
    parms <- unique(d$parameter)
    tmp <- vector("list", length(parms))
    for(i in seq_along(parms)) {
        para <- parms[[i]]
        ind <- d$parameter == para
        se <- d[ind, "se"]
        estimate <- d[ind, "estimate"]
        pval <- d[ind, "pval"]
        df <- d[ind, "df"]
        df_bw <- d[ind, "df_bw"]

        if(any(!is.na(pval))) {
            Satt_NA <- mean(is.na(df))
        } else Satt_NA <- NA

       # para <- i
        #theta_i <- theta[[i]]
        pvals_bw <- get_p_val_df(t = estimate/se,
                                 df = df_bw,
                                 parameter = para)
        theta_i <- theta[[i]]
        tmp[[i]] <- data.frame(
                    parameter = para,
                    M_est = mean(estimate),
                    theta = theta_i,
                    M_se = mean(se),
                    SD_est = sd(estimate),
                    Power = mean(get_cover(estimate, se, alpha = alpha)),
                    Power_bw = mean(pvals_bw < alpha),
                    Power_satt = mean(pval < alpha, na.rm = TRUE),
                    Satt_NA = Satt_NA
                )


    }
    FE <- do.call(rbind, tmp)

    FE
}
summarize_CI <- function(res, theta) {
    d <- res$FE
    parms <- unique(res$FE$parameter)
    tmp <- vector("list", length(parms))
    for(i in seq_along(parms)) {
        para <- parms[[i]]
        ind <- d$parameter == para
        CI_lwr <- d[ind, "CI_lwr"]
        CI_upr <- d[ind, "CI_upr"]
        CI_wald_lwr <- d[ind, "CI_wald_lwr"]
        CI_wald_upr <- d[ind, "CI_wald_upr"]

        tmp[[i]] <- data.frame(
                    parameter = para,
                    CI_cover = mean(CI_lwr < theta &
                                         CI_upr > theta, na.rm=TRUE),
                    CI_Wald_cover = mean(CI_wald_lwr < theta &
                                             CI_wald_upr > theta, na.rm=TRUE),
                    CI_NA = mean(is.na(CI_lwr)))
    }
    CI_cov <- do.call(rbind, tmp)

    CI_cov
}
summary_.plcp_sim  <- function(res, paras, alpha) {
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

    RE <- summarize_RE(res, theta = RE_params)

    # falseconv
    false_conv <-
        res$RE[res$RE$parameter == "cluster_slope", "vcov"]
    false_conv <- mean(is_approx(false_conv, 0))


    theta = list(
        "(Intercept)" = paras$fixed_intercept,
        "treatment" = 0,
        "time" = paras$fixed_slope,
        "time:treatment" = get_slope_diff(paras) / paras$T_end
    )
    FE <- summarize_FE(res = res,
                       theta = theta,
                       alpha = alpha)

    if(all(is.na(res$FE$pval))) {
        Satt_NA <- 0
    } else {
        Satt_NA <- FE[FE$parameter == "time:treatment", "Satt_NA"]
    }

    FE$Satt_NA <- NULL


    CI_NA <- 0
    if ("CI_lwr" %in% colnames(res$FE)) {
        CI_cov <- summarize_CI(res, theta[["time:treatment"]])
        FE$CI_Cover <- CI_cov$CI_cover
        CI_NA <- CI_cov[CI_cov$parameter == "time:treatment", "CI_NA"]
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

    convergence <- mean(res$convergence)

    list("RE" = RE,
         "FE" = FE,
         "tot_n" = res$tot_n,
         "convergence" = convergence,
         "Satt_NA" = Satt_NA,
         "CI_NA" = CI_NA)
}


# multi-sim ---------------------------------------------------------------


# as.plcp_sim <- function(object) {
#     UseMethod("as.plcp_sim")
# }
#
# as.plcp_sim.data.frame <- function(object) {
#     object <- as.list(object)
#     class(object) <- append("plcp", class(object))
#
#     object
# }


#' Summarize simulations based on a combination of multiple parameter values
#'
#' @param object A multiple simulation object created with
#' \code{\link{simulate.plcp_multi}}
#' @param para The name of the fixed or random effect that should be summarized.
#' @param type Specifies what type of effect \code{para} is; can be either
#' "random" or "fixed".
#' @param model Specifies which model should be summarized. Used when models are
#' simulated from formulas named "correct" and "wrong".
#' @param alpha Indicates the significance level. Default is 0.05 (two-tailed),
#' one-tailed tests are not yet implemented.
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
                                   model = "correct",
                                   alpha = 0.05,
                                   ...) {
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
        out <- lapply(res, summary_fixed.plcp_multi_sim, para, model, alpha = alpha)
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


#' Print method for \code{summary.plcp_multi_sim}-objects
#' @param x An object of class \code{plcp_multi_sim_summary}
#' @param digits number of significant digits.
#' @param ... Optional arguments.
#' @method print plcp_multi_sim_summary
#' @export
print.plcp_multi_sim_summary <- function(x, digits = 2, ...) {
    model <- attr(x, "model")
    type <- attr(x, "type")
    nsim <- attr(x, "nsim")
    x <- as.data.frame(x)
    cat("Model: ", model, "| Type:", type, "\n")
    cat("---\n")
    print(x, digits = digits)
    cat("---\n")
    cat("nsim: ", nsim, "\n")
    if(any(x$is_NA > 0)) warning("Some simulations had NA estimates that was removed.", call. = FALSE)
    invisible(x)
}

## random effect
summary_fixed.plcp_multi_sim <- function(res, para, model, alpha) {

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

    power <- mean(get_cover(out$estimate,
                            out$se,
                            alpha = alpha))
    power_bw <- get_p_val_df(t = out$estimate/out$se,
                             df = out$df_bw,
                             parameter = para)
    Satt_NA <-  mean(is.na(out$df))

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
            Power_bw = mean(power_bw < alpha),
            Power_satt = mean(out$pval < alpha, na.rm = TRUE),
            Satt_NA = Satt_NA

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
        "subject_intercept" = res$paras$sigma_subject_intercept^2,
        "subject_slope" = res$paras$sigma_subject_slope^2,
        "subject_intercept" = res$paras$sigma_subject_intercept^2,
        "cor_subject" = res$paras$cor_subject,
        "cluster_intercept" = res$paras$sigma_cluster_intercept^2,
        "cluster_slope" = res$paras$sigma_cluster_slope^2,
        "cor_cluster" = res$paras$cor_cluster,
        "error" = res$paras$sigma_error^2
    )

    out <- res$res[[model]]$RE
    out <- out[out$parameter == para,]

    vcov <- out$vcov
    est <- mean(vcov, na.rm = TRUE)
    out <- data.frame(
        parameter = para,
        M_est = est,
        theta = theta,
        est_rel_bias = get_RB(est, theta),
        prop_zero = mean(is_approx(vcov, 0), na.rm=TRUE),
        is_NA = mean(is.na(vcov))
    )
    if (para %in% c("cor_subject", "cor_cluster")) {
        out$prop_zero <- mean(is_approx(abs(vcov), 1), na.rm=TRUE)
    }

    out

}


get_RB <- function(est, theta) {
    if(theta == 0) {
        RB <- (est - theta)
    } else {
        RB <- (est - theta) / theta
    }
    RB
}

# power -------------------------------------------------------------------
get_cover <- function(est, se, alpha) {
    # CI cover
    abs(est) - qnorm(1 - (alpha/2)) * se > 0
}

get_p_val_df <- function(t, df, parameter) {
    if (unique(parameter) %in% c("time:treatment", "time")) {
        p <- (1 - pt(abs(t), df)) * 2
    } else
        p <- NA

    p
}
