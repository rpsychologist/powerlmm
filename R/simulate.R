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



# sim formula -------------------------------------------------------------

#' Create a simulation formula
#'
#' @param formula A \code{character} containing a \pkg{lme4} formula.
#' @param data_transform Optional; a \code{function} that applies a transformation
#' to the data during each simulation.
#' @param test A \code{character} vector indicating which parameters should be tested.
#' Only applies to tests using Satterthwaite \emph{dfs}, or when calculating confidence intervals.
#'
#' @details
#'
#' It is possible to fit model without any random effects. If no random effects is specified
#' the model is fit using \code{lm()}.
#'
#' @return Object with class \code{plcp_sim_formula}
#' @seealso \code{\link{sim_formula_compare}}, \code{\link{transform_to_posttest}}
#' @export
#'
#' @examples
#' # 2-lvl model
#' f <- sim_formula("y ~ treatment * time + (1 + time | subject)")
#'
#' # ANCOVA using 'data_transform'
#' f <- sim_formula("y ~ treatment + pretest",
#'                  data_transform = transform_to_posttest,
#'                  test = "treatment")
#'
sim_formula <- function(formula, data_transform = NULL,  test = "time:treatment", ...) {
    UseMethod("sim_formula")
}
#' @export
sim_formula.default <- function(formula, data_transform = NULL, test = "time:treatment", ...) {

     x <- list("formula" = formula,
              "data_transform" = data_transform,
              "data_transform_lab" = substitute(data_transform),
              "test" = test, ...)

    class(x) <- append(class(x), "plcp_sim_formula")
    x
}
#' Compare multiple simulation formulas
#'
#' This functions allows comparing multiple models fit to the same data set
#' during simulation.
#'
#' @param ... Named formulas that should be compared, see \emph{Examples}.
#'
#' @return Object with class \code{plcp_compare_sim_formula}
#' @seealso \code{\link{sim_formula}}
#' @export
#'
#' @examples
#'
#' # Formulas can be a named character
#' # uses the defaults 'sim_formula()'
#' f <- sim_formula_compare("m0" = "y ~ time * treatment + (1 | subject)",
#'                          "m1" = "y ~ time * treatment + (1 + time | subject)")
#'
#' # Can also use sim_formula()
#' f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
#' f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
#' f <- sim_formula_compare("m0" = f0, "m1" = f1)
#'
sim_formula_compare <- function(...) {
    x <- list(...)
    if(any(is.null(names(x)))) stop("formula(s) must have a name.", call. = FALSE)
    for(i in seq_along(x)) {
        f <- x[[i]]
        if(is.character(f)) x[[i]] <- sim_formula(f)
    }

    class(x) <- append(class(x), "plcp_compare_sim_formula")
    x
}


#' Print method for simulation formulas
#'
#' @param x A formula object.
#' @param ... Not used
#'
#' @export
#'
#' @method print plcp_sim_formula
print.plcp_sim_formula <- function(x, ...) {
    cat("# Simulation formula\n")
    f <- x$formula
    data_transform <- x$data_transform_lab
    test <- paste(x$test, collapse = "', '")
    if(!is.null(data_transform))
        data_transform <- paste0("\n  data_transform: '", data_transform, "'")
    cat(paste("         formula: '", f,
              "'\n            test: '", test, "'",
              data_transform,
              sep = ""), sep = "\n")
    cat("\n")

}

#' @rdname print.plcp_sim_formula
#' @method print plcp_compare_sim_formula
#' @export
print.plcp_compare_sim_formula <- function(x, ...) {

    cat("# Simulation formulas\n")
    lapply(seq_along(x), .print_plcp_sim_formula, x = x)

}
.print_plcp_sim_formula <- function(i, x) {
    lab <- names(x)[i]
    f <- x[[i]]$formula
    data_transform <- x[[i]]$data_transform_lab
    test <- paste(x[[i]]$test, collapse = "', '")
    if(!is.null(data_transform))
        data_transform <- paste0("\n             data_transform: '", data_transform, "'")
    cat(paste("#", i, "    '", lab,"': '", f,
              "'\n             test: '", test, "'",
              data_transform,
              sep = ""), sep = "\n")
    cat("\n")
}

# data transform helpers --------------------------------------------------


#' Helper to transform the simulated longitudinal \code{data.frame}
#'
#' This is en example of a data transformation applied during simulation.
#' It takes the longitudinal data and transforms it into a pretest-posttest
#' model in wide format. Useful if you want to compare the longitudinal LMM with
#' e.g. AN(C)OVA models.
#'
#' @param data a \code{data.frame} created using \code{\link{simulate_data}}
#'
#' @return a \code{data.frame} with \code{y} now only includes the posttest values.
#' Also includes three new columns:
#' \itemize{
#'    \item \code{pre} subject-level pretest scores.
#'    \item \code{pre_cluster} cluster-level pretest scores.
#'    \item \code{pre_subject_c} subject-level pretest scores center
#'    around the cluster-level pretest.
#' }
#' @export
#'
#' @seealso \code{\link{simulate.plcp}}, \code{\link{study_parameters}}
#'
#' @examples
#'
#' # Compare longitudinal 3-level model to 2-level model
#' # fit to just the posttest data
#' #
#' # Both models are fit to the same dataset during simulation.
#' p <- study_parameters(n1 = 11,
#'                       n2 = 20,
#'                       n3 = 3,
#'                       icc_pre_subject = 0.5,
#'                       icc_pre_cluster = 0.1,
#'                       icc_slope = 0.05,
#'                       var_ratio = 0.03)
#'
#' # simulation formulas
#' # analyze as a posttest only 2-level model
#' f_pt <- sim_formula("y ~ treatment + (1 | cluster)",
#'                  test = "treatment",
#'                  data_transform = transform_to_posttest)
#'
#' # analyze as 3-level longitudinal
#' f_lt <- sim_formula("y ~ time*treatment +
#'                          (1 + time | subject) +
#'                          (1 + time | cluster)")
#'
#' f <- sim_formula_compare("posttest" = f_pt,
#'                          "longitudinal" = f_lt)
#' \dontrun{
#' res <- simulate(p,
#'                 formula = f,
#'                 nsim = 2000,
#'                 cores = parallel::detectCores(),
#'                 satterthwaite = TRUE)
#' summary(res)
#'}
#'
transform_to_posttest <- function(data) {
    tmp <- data[data$time == max(data$time), ]
    tmp$pretest <- data[data$time == min(data$time), "y"]
    tmp$pretest_cluster <- tapply(tmp$pretest, tmp$cluster, mean)[tmp$cluster]
    tmp$pretest_subject_c <- tmp$pretest - tmp$pretest_cluster

    tmp
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
simulate_data.plcp_nested <- function(paras, n = NULL) {

    # get sim data func
    if(paras$custom_model) {
        message("custom")
    } else {
        .sim_data_func <- .simulate_3lvl_data
    }


    if (is.data.frame(paras))
        paras <- as.list(paras)
    if(is.null(paras$prepared)) {
        tmp <- prepare_paras(paras)
    } else tmp <- paras
    paras <- tmp$control
    paras_tx <- tmp$treatment

    slope_diff <- get_slope_diff(paras) / paras$T_end
    paras_tx$fixed_slope <- paras_tx$fixed_slope + slope_diff

    # replace NA
    paras[is.na(paras)] <- 0
    paras_tx[is.na(paras_tx)] <- 0

    d_tx <- do.call(.sim_data_func, paras_tx)
    d_c <- do.call(.sim_data_func, paras)

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
simulate_data.plcp <- function(paras, ...) {
    NextMethod("simulate_data")
}

# Crossed data simulate  --------------------------------------------------
#' @rdname simulate_data
#' @export
simulate_data.plcp_crossed <- function(paras, n = NULL) {
   # if (is.data.frame(paras))
    #    paras <- as.list(paras)
    if(is.null(paras$prepared)) {
        tmp <- prepare_paras(paras)
    } else {
        tmp <- paras
        paras <- paras$control
    }


    # doesn't matter if tmp$control or tx
    slope_diff <- get_slope_diff(tmp$control) / paras$T_end
    paras$effect_size <- NULL

    paras$n2 <- list(control = tmp$control$n2,
                     treatment = tmp$treatment$n2)
    paras$fixed_slope_time_tx <- slope_diff

    # replace NA
    paras[is.na(paras)] <- 0

    d <- simulate_3lvl_data_crossed(paras)

    # # drop outs

    if (is.list(paras$dropout) |
        is.function(paras$dropout) |
        length(paras$dropout) > 1) {

        miss_c <- create_dropout_indicator(tmp$control)
        miss_tx <- create_dropout_indicator(tmp$treatment)

        d <- add_NA_values_from_indicator(d, c(miss_tx, miss_c))
    }

    d
}


#' Perform a simulation study using a \code{study_parameters}-object
#'
#' @param object An object created by \code{\link{study_parameters}}.
#' @param nsim The number of simulations to run.
#' @param seed Currently ignored.
#' @param formula Model formula(s) used to analyze the data, see \emph{Details}.
#' Should be created using \code{\link{sim_formula}}. It is also possible to compare multiple
#' models, e.g. a correct and a misspecified model, by combining the formulas using \code{\link{sim_formula_compare}}.
#' See \emph{Examples}. If \code{NULL} the formula is made automatically,
#' using \code{\link{create_lmer_formula}}, which does not support objects with multiple
#' simulation setups.
#' @param satterthwaite Logical; if \code{TRUE} Satterthwaite's degrees of freedom
#' approximation will be used when computing \emph{p}-values. This is implemented using
#' the \code{lmerTest}-package. See \emph{Details}.
#' @param CI Logical; if \code{TRUE} coverage rates for 95 \% confidence intervals
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
#' @details
#'
#' See also \code{vignette("simulations", package = "powerlmm")} for a tutorial.
#'
#' \strong{Model formula}
#'
#' If no data transformation is used, the available model terms are:
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
#' See \emph{Examples} and the simulation-vignette for formula examples. For
#' \code{object}s that contain a single study setup, then the lmer formula
#' can be created automatically using \code{\link{create_lmer_formula}}.
#'
#' \strong{Satterthwaite's approximation, and CI coverage}
#'
#' To decrease the simulation time the default is to only calculate Satterthwaite's \emph{dfs}
#' and the CIs' coverage rates for the test of 'time:treatment'-interaction. This can be
#' changed using the argument \code{test} in \code{\link{sim_formula}}.
#'
#' Confidence intervals are both calculated using profile likelihood and by
#' the Wald approximation, using a 95 \% confidence level.
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
#' @seealso \code{\link{sim_formula}}, \code{\link{sim_formula_compare}}, \code{\link{summary.plcp_sim}}, \code{\link{simulate_data}}
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
#' f <- sim_formula("y ~ treatment * time + (1 + time | subject)")
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
#' # Three-level (nested) ------------------------------------------------------
#' p <- study_parameters(n1 = 10,
#'                       n2 = 20,
#'                       n3 = 4,
#'                       sigma_subject_intercept = 1.44,
#'                       icc_pre_cluster = 0,
#'                       sigma_subject_slope = 0.2,
#'                       icc_slope = 0.05,
#'                       sigma_error = 1.44,
#'                       effect_size = 0)
#'
#' ## compare correct and miss-specified model
#' f0 <- "y ~ treatment * time + (1 + time | subject)"
#' f1 <- "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"
#' f <- sim_formula_compare("correct" = f1,
#'                          "wrong" = f0)
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
#' ## Compare random effects using LRT,
#' ## summarise based on best model from each sim
#' summary(res,
#'         model_selection = "FW",
#'         LRT_alpha = 0.1,
#'         para = "treatment:time")
#'
#' # Partially nested design ---------------------------------------------------
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
#' f <- sim_formula("y ~ treatment * time + (1 + time | subject) +
#'                   (0 + treatment:time | cluster)")
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
#' # Run multiple designs  -----------------------------------------------------
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
#' f0 <- "y ~ treatment * time + (1 + time | subject)"
#' f1 <- "y ~ treatment * time + (1 + time | subject) + (0 + time | cluster)"
#' f <- sim_formula_compare("correct" = f1,
#'                          "wrong" = f0)
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
#' summary(res, para = "time:treatment", model = "correct")
#'
#' # Summarize cluster-level random slope  for n3 = c(2, 4, 6) for 'correct' model
#' summary(res, para = "cluster_slope", model = "correct")
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

    if(is.null(formula)) formula <- sim_formula(create_lmer_formula(object))
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
        if(inherits(formula, "plcp_compare_sim_formula")) {
            class(out) <- append(class(out), "plcp_sim_formula_compare")
        }
         else {
            class(out) <- append(class(out), "plcp_sim")
        }
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
        attr(res, "paras") <- object
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
                           tot_n = tot_n,
                           sim = sim)

    res
}


# Checks ------------------------------------------------------------------

# checks


check_formula <- function(formula) {
    if(!inherits(formula, "plcp_sim_formula") & !inherits(formula, "plcp_compare_sim_formula")) {
        stop("`formula` should be created using `sim_formula` or `sim_formula_compare`", call. = FALSE)
    }
    #if (!all(names(formula) %in% c("correct", "wrong")))
    #    stop("Formula names must be either 'correct' or 'wrong'")

    if(inherits(formula, "plcp_sim_formula")) formula <- list("default" = formula)
    if (inherits(formula, "plcp_compare_sim_formula") & length(names(formula)) > 1) {
        if (length(names(formula)) != length(unique(names(formula)))) {
            stop(
                "Formulas should have unique names."
            )
        }
    }

    #for (f in formula)
    #    check_formula_terms(f)

    formula
}
check_formula_terms <- function(f) {
    f <- as.formula(f$formula)

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

fit_model <- function(formula, data, ...) {
    UseMethod("fit_model")
}
fit_model.default <- function(formula, data, ...) {
    # LMM or OLS

    formula <- as.formula(formula$formula)
    if(is.null(lme4::findbars(formula))) {
        fit <- tryCatch(
            #do.call(lme4::lmer, list(formula=f, data=d))
            fit <- stats::lm(formula = formula, data = data)
        )
    } else {
        fit <- tryCatch(
            #do.call(lme4::lmer, list(formula=f, data=d))
            fit <- lme4::lmer(formula = formula, data = data)
        )
    }
}

analyze_data <- function(formula, d) {
    fit <-
        lapply(formula, function(f) {
            #if(inherits(f, "plcp_sim_formula")) {
                if(is.function(f$data_transform))
                    d <- f$data_transform(d)

            #}
            fit <- fit_model(f,
                             data = d)

           list("fit" = fit,
                "test" = f$test,
                "formula" = f)
        })

    fit
}

extract_results <- function(fit, CI = FALSE, satterthwaite = FALSE, df_bw, tot_n, sim) {
    lapply(fit, extract_results_,
           CI = CI,
           satterthwaite = satterthwaite,
           df_bw = df_bw,
           tot_n = tot_n,
           sim = sim)
}
extract_random_effects <- function(fit) {
    UseMethod("extract_random_effects")
}
extract_random_effects.lm <- function(fit) {
    data.frame(grp = "Residual",
               var1 = NA,
               var2 = NA,
               vcov = stats::sigma(fit)^2,
               parameter = "Residual_NA",
               stringsAsFactors = FALSE)
}
extract_random_effects.lmerMod <- function(fit) {
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

#' @export
add_p_value <- function(fit, test, ...) {
    UseMethod("add_p_value")
}
#' @importFrom utils packageVersion
add_p_value.lmerMod <- function(fit, test, satterthwaite, df_bw = NULL, ...) {
    if(satterthwaite) {
        tmp <- tryCatch(summary(fit))
        tmp <- tmp$coefficients

        ff <- rownames(tmp)
        # satterth dfs only for time:treatment

        satterth_term <- any(ff %in% test)
        ## satterthwaite
        if(satterth_term) {
            ind <- which(ff %in% test)

            L <- rep(0, length(ff))
            L <- rep(list(L), length(ind))

            for(i in seq_along(L)) {
                L[[i]][ind[i]] <- 1
            }

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
            if(any(is.na(p[ind]))) {
                ind <- which(is.na(p[ind]))
                tval <- tmp[ind, "t value"]
                p[ind] <- 2*(1 - pt(abs(tval), df = df_bw))
            }

            res <- list("df" = df,
                        "p" = p)
        }
    } else {
        res <- list("df" = NA,
                    "p" = NA)
    }


    res
}
add_p_value.lm <- function(fit, test, ...) {

    tmp <- summary(fit)
    tmp <- as.data.frame(tmp$coefficients)
    ff <- rownames(tmp)
    ind <- which(ff %in% test)
    # satterth dfs only for time:treatment

    if(any(ff %in% test)) {
        df <- rep(NA, length(ff))
        p <- rep(NA, length(ff))
        df[ind] <- fit$df.residual
        p[ind] <- tmp$`Pr(>|t|)`[ind]

        list("df" = df,
             "p" = p)
    } else {
        list("df" = NA,
             "p" = NA)
    }

}
fix_sath_NA_pval <- function(x, df) {
    ind <- is.na(x$pval)
    tmp <- x[ind, ]
    t <- with(tmp, estimate/se)
    pval <- 2*(1 - pt(abs(t), df = df))

    x[ind, "pval"] <- pval

    x
}

#' @export
get_fixef <- function(fit, ...) {
    UseMethod("get_fixef")
}
get_fixef.default <- function(fit, test, df_bw, satterthwaite, ...) {

    # need to know in which order time and treatment was entered
    # then adjust 'fit$test' to match
    FE_coefs <- get_fixef_coef(fit, ...)
    rnames <- names(FE_coefs)
    TbT <- c("time:treatment", "treatment:time")
    ind <- TbT %in% rnames
    if(any(ind) & any(test %in% TbT)) {
        TbT_model <- TbT[TbT %in% rnames]
        test[test %in% TbT] <- TbT_model
    }

    if(any(!test %in% rnames)) stop("At least of the values in 'test' do no match the model's parameters: ", paste(rnames, collapse = ", "), call. = FALSE)

    se <- sqrt(diag(vcov(fit)))
    tmp_p <- add_p_value(fit = fit,
                         test = test,
                         satterthwaite = satterthwaite,
                         df_bw = df_bw)
    FE <- data.frame(
        "estimate" = FE_coefs,
        "se" = se,
        "pval" = tmp_p$p,
        "df" = tmp_p$df
    )

    rownames(FE) <- NULL
    FE <- cbind(data.frame(parameter = rnames,
                           stringsAsFactors = FALSE),
                FE)

    FE[FE$parameter %in% test, "df_bw"] <- df_bw

    FE
}

#' @export
get_fixef_coef <- function(fit, ...) {
    UseMethod("get_fixef_coef")
}

get_fixef_coef.lmerMod <- function(fit, ...) {
    lme4::fixef(fit)
}
get_fixef_coef.lm <- function(fit, ...) {
    stats::coef(fit)
}

#' @export
get_convergence <- function(fit) {
    UseMethod("get_convergence")
}
get_convergence.lmerMod <- function(fit) {
    is.null(fit@optinfo$conv$lme4$code)
}
get_convergence.lm <- function(fit) {
    TRUE
}

#' @export
get_CI <- function(fit, test, FE, ...) {
    UseMethod("get_CI")
}
get_CI.default <- function(fit, test, FE, ...) {
    CI <- tryCatch(confint(fit, parm = test),
                   error = function(e) NA)
    CI_wald <-
        tryCatch(confint(fit, method = "Wald", parm = test),
                 error = function(e) NA)
    CIs <- test
    if(all(is.na(CI))) {
        FE[FE$parameter %in% CIs, "CI_lwr"] <- NA
        FE[FE$parameter %in% CIs, "CI_upr"] <- NA
    } else {
        FE[FE$parameter %in% CIs, "CI_lwr"] <- CI[CIs, 1]
        FE[FE$parameter %in% CIs, "CI_upr"] <- CI[CIs, 2]
    }
    if(all(is.na(CI_wald))) {
        FE[FE$parameter %in% CIs, "CI_wald_lwr"] <- NA
        FE[FE$parameter %in% CIs, "CI_wald_upr"] <- NA
    } else {
        FE[FE$parameter %in% CIs, "CI_wald_lwr"] <- CI_wald[CIs, 1]
        FE[FE$parameter %in% CIs, "CI_wald_upr"] <- CI_wald[CIs, 2]
    }
    FE
}

#' @export
get_LL <- function(fit) {
    UseMethod("get_LL")
}
get_LL.default <- function(fit) {
    ll <- stats::logLik(fit, REML = TRUE)
    df <- attr(ll, "df")
    list(ll, df)
}

extract_results_ <- function(fit, CI, satterthwaite,  df_bw, tot_n, sim) {

    FE <- get_fixef(fit = fit$fit,
                    test = fit$test,
                    satterthwaite = satterthwaite,
                    df_bw = df_bw,
                    formula = fit$formula)

    if (CI) {
        FE <- get_CI(fit = fit$fit,
                     test = fit$test,
                     FE = FE)
    }
    RE <- extract_random_effects(fit$fit)

    conv <- get_convergence(fit$fit)

    # save for postprocess LRT test
    ll <- get_LL(fit$fit)
    RE$sim <- sim
    FE$sim <- sim

    out <- list("RE" = RE,
         "FE" = FE,
         "logLik" = as.numeric(ll$ll),
         "df" = ll$df,
         "tot_n" = tot_n,
         "conv" = conv)
    class(out) <- append(class(out), class(fit$fit))

    out
}

rename_rr_ <- function(.x, match, new) {
    .x[.x$parameter %in% match, "parameter"] <- new

    .x
}

#' @export
rename_random_effects <- function(.x, crossed = FALSE) {
    UseMethod("rename_random_effects")
}
#' @export
rename_random_effects.default <- function(.x, crossed = FALSE) {
    if(crossed) {
        .x <- rename_rr_(
            .x,
            match = c(
                "cluster_time",
                "cluster.1_time"
            ),
            new = "cluster_slope"
        )
        .x <- rename_rr_(
            .x,
            match = c(
                "cluster_time:treatment",
                "cluster_treatment:time",
                "cluster.1_treatment:time",
                "cluster.1_time:treatment"
            ),
            new = "cluster_slope_tx"
        )

    } else {
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
    }

    .x <- rename_rr_(.x,
                     match = "Residual_NA",
                     new = "error")

    .x <- rename_rr_(.x,
                     match = c("subject_(Intercept)"),
                     new = "subject_intercept")

    if(crossed) {
        .x <- rename_rr_(.x,
                         match = c("cluster_(Intercept)"),
                         new = "cluster_intercept")
        .x <- rename_rr_(.x,
                         match = c("cluster_treatment"),
                         new = "cluster_intercept_tx")
    } else {
        .x <- rename_rr_(.x,
                         match = c("cluster_(Intercept)",
                                   "cluster_treatment"),
                         new = "cluster_intercept")
    }



    .x <- rename_rr_(.x,
                     match = c("subject.1_time", "subject_time"),
                     new = "subject_slope")


    .x <- rename_rr_(.x,
                     match = "subject_(Intercept)_time",
                     new = "cor_subject")

    if(crossed) {
        .x <- rename_rr_(.x,
                         match = c("cluster_(Intercept)_time"),
                         new = "cor_cluster_intercept_slope")
        .x <- rename_rr_(.x,
                         match = c("cluster_(Intercept)_treatment:time",
                                   "cluster_(Intercept)_time:treatment",
                                   "cluster_treatment_treatment:time"),
                         new = "cor_cluster_intercept_slope_tx")
    } else {
        .x <- rename_rr_(.x,
                         match = c("cluster_(Intercept)_time",
                                   "cluster_(Intercept)_treatment:time",
                                   "cluster_(Intercept)_time:treatment",
                                   "cluster_treatment_treatment:time"),
                         new = "cor_cluster")
    }


    # crossed effects
    ## cor intercept, time
    .x <- rename_rr_(.x,
                     match = "cluster_time_treatment",
                     new = "cor_cluster_slope_intercept_tx")
    ## cor intercept, treatment
    .x <- rename_rr_(.x,
                     match = "cluster_(Intercept)_treatment",
                     new = "cor_cluster_intercept_intercept_tx")

    ## cor tx, time:tx
    .x <- rename_rr_(.x,
                     match = "cluster_treatment_time:treatment",
                     new = "cor_cluster_intercept_tx_slope_tx")

    ## cor time, time:tx
    .x <- rename_rr_(.x,
                     match = "cluster_time_time:treatment",
                     new = "cor_cluster_slope_slope_tx")


}
order_model_lists <- function(models, RE, FE, ll, df, tot_n, convergence) {
    x <- list()
    for (i in seq_along(models)) {
        x[[i]] <- list("RE" = RE[[i]],
                       "FE" = FE[[i]],
                       "logLik" = ll[[i]],
                       "df" = df[[i]],
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
    crossed <- inherits(res$paras, "plcp_crossed")


    for(i in seq_along(RE)) {
        tmp <- RE[[i]]
        class(tmp) <- append(class(tmp), grep("plcp", class(res$formula[[i]]), value = TRUE))
        RE[[i]] <- tmp
    }
    RE <- lapply(RE, rename_random_effects, crossed = crossed)
    ll <- lapply(models, munge_results_, x, "logLik")
    df <- lapply(models, munge_results_, x, "df")

    x <- order_model_lists(models,
                           RE = RE,
                           FE = FE,
                           ll = ll,
                           df = df,
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

.print_effect <- function(effect, lab, digits) {
    if(nrow(effect) >= 1) {
        if(c("model" %in% colnames(effect))) {
            fix_eff <- paste(as.character(unique(effect$parameter)), collapse = "', '")
            tmp <- effect[ , colnames(effect) %in% c("model"), drop = FALSE]
            cat("\n", lab, ": '", fix_eff, "'\n\n", sep = "")
        } else {
            cat("\n", lab, " \n\n", sep = "")
            tmp <- effect[ , colnames(effect) %in% c("parameter", "model"), drop = FALSE]
        }
        effect <- effect[, !colnames(effect) %in% c("parameter", "model"), drop = FALSE]
        effect <- signif(effect, digits)
        effect <- cbind(tmp, effect)
        effect[t(do.call(rbind, lapply(effect, is.nan)))] <- "."
        effect[t(do.call(rbind, lapply(effect, is.na)))] <- "."
        print.data.frame(effect,
                         row.names = FALSE,
                         digits = digits,
                         quote = FALSE)
    }
}

print_model <- function(i, x, digits = 2) {
    models <- names(x)
    cat("Model: ", models[i], "\n")
    .print_effect(x[[i]]$RE,
                  lab = "Random effects",
                  digits = digits)
    .print_effect(x[[i]]$FE,
                  lab = "Fixed effects",
                  digits = digits)
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
    if(!is.na(Satt_NA) && Satt_NA > 0) {
        message("[Model: ", mod_name, "] ", round(Satt_NA, 4) * 100, "% of the Satterthwaite calculations failed")
    }
    if(!is.na(CI_NA) && CI_NA > 0) {
        message("[Model: ", mod_name, "] ", round(CI_NA, 4) * 100, "% of the profile likelihood CIs failed")
    }
    if(!is.na(non_convergence) && non_convergence > 0) {
        message("[Model: ", mod_name, "] ", round(non_convergence, 4) * 100, "% of the models threw convergence warnings")
    }
}

#' Print method for \code{summary.plcp_sim}-objects
#' @param x An object of class \code{plcp_sim_summary}
#' @param verbose \code{logical}; indicates if additional information
#' should be printed (default is \code{TRUE}).
#' @param digits number of significant digits.
#' @param ... Optional arguments.
#' @method print plcp_sim_summary
#' @export
print.plcp_sim_summary <- function(x, verbose = TRUE, digits = 2, ...) {
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
    if(verbose) {
        cat("Number of simulations:", res$nsim, " | alpha: ", res$alpha)
        cat("\nTime points (n1): ", res$paras$n1)
        cat(n2_lab, n2)
        cat("\nTotal number of subjects: ", tot_n, "\n")
        lapply(seq_along(x), print_test_NA, x = x)
        if("model_selected" %in% names(res)) {
            cat("---\nResults based on LRT model comparisons, using direction: ",
                res$model_direction, " (alpha = ", res$LRT_alpha, ")\n", sep = "")
            print(res$model_selected)
        }
        data_bool <- vapply(res$formula, function(x) is.function(x$data_transform), logical(1))
        if(any(data_bool)) cat("---\nAt least one of the models applied a data transformation during simulation,\n",
                               "summaries that depend on the true parameter values will no longer be correct,\n",
                               "see 'help(summary.plcp_sim)'\n", sep = "")
    }
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

    ff <- x[[1]]$formula
    if(inherits(ff, "plcp_compare_sim_formula")) {
        print(ff)
    } else if(inherits(ff[[1]], "plcp_sim_formula")) {
        print(ff[[1]])
    }

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
    print(x$formula[[1]])
    cat("# Object size:", format(object.size(x), units = "auto"))
    cat("\n")
    invisible(x)
}

#' @rdname print.plcp_sim
#' @method print plcp_sim_formula_compare
#' @export
print.plcp_sim_formula_compare <- function(x, ...) {
    cat(
        "# A 'plcp_sim'-object containing",
        x$nsim,
        "simulations.",
        "Use summary() to view results."
    )
    cat("\n")
    print(x$formula)
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
#' @param para Selects a parameter to return. Default is \code{NULL},
#' which returns all parameters. If multiple model formulas are compared a named list can be
#' used to specify different parameters per model.
#' @param ... Currently not used
#'
#' @details
#'
#' \strong{Model selection}
#'
#' It is possible to summarize the performance of a data driven model selection strategy
#' based on the formulas used in the simulation (see \code{\link{sim_formula_compare}}).
#' The two model selection strategies are:
#' \itemize{
#'   \item \code{FW}: Forward selection of the models. Starts with the first model formula and
#'   compares it with the next formula. Continues until the test of M_i vs M_{i + 1} is non-significant,
#'   and then picks M_i. Thus if three models are compared, and the comparison of M_1 vs M_2 is non-significant, M_3
#'   will not be tested and M_1 is the winning model.
#'  \item \code{BW}: Backward selection of the models. Starts with the last model formula and
#'   compares it with the previous formula. Continues until the test of M_i vs M_{i - 1} is significant or until
#'   all adjacent formulas have been compared. Thus if three models are compared, and the comparison of M_3 vs M_2 is non-significant,
#'   M2 vs M1 will be tested and M2 will be picked if significant, and M1 if not.
#' }
#'
#' The model comparison is performed using a likelihood ratio test based the REML criterion. Hence, it assumed you are comparing models
#' with the same fixed effects, and that one of the models is a reduced version of the other (nested models). The LRT test is done as a
#' post-processing step, so \code{model_selection} option will not re-run the simulation. This also means that different alpha levels
#' for the LRTs can be investigated without re-running the simulation.
#'
#' \strong{Data transformation}
#'
#' If the data has been transformed \code{sim_formula(data_transform = ...)}, then
#' true parameter values (\code{theta}s shown in the summary will most likely no longer
#' apply. Hence, relative bias and CI coverage will be in relation to the original model.
#' However, the empirical estimates will be summarized correctly, enabling investigation of
#' power and Type I errors using arbitrary transformations.
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
summary.plcp_sim <- function(object, model = NULL, alpha = 0.05, para = NULL, ...) {
    res <- object
    x <- lapply(res$res, summary_.plcp_sim,
                paras = res$paras,
                alpha = alpha,
                ...)
    x <- list(summary = x,
              nsim = res$nsim,
              paras = res$paras,
              tot_n = x[[1]]$tot_n,
              alpha = alpha,
              formula = res$formula)

    if(!is.null(para)) {
        #check_para <- vapply(x$summary, function(x) para %in% x$FE$parameter, logical(1))
        #if(!all(check_para)) stop("The parameter: '", para, "' does not exist in all models.", call. = FALSE)
        nr <- length(x$summary)
        FE <- vector("list", nr)
        RE <- vector("list", nr)
        for(i in 1:nr) {
            mod <- names(x$summary)[i]
            tmp <- x$summary[[i]]
            if(is.list(para)) {
                pp <- para[[mod]]
            } else if(is.character(para)) {
                if(length(para) > 1) stop("'para' must have length equal to 1. ",
                                          "If you want different parameters per model ",
                                          "formula use a named list.", call. = FALSE)
                pp <- para
            }
            FE[[i]] <- tmp$FE[tmp$FE$parameter == pp, ]
            RE[[i]] <- tmp$RE[tmp$RE$parameter == pp, ]
            if(is.list(para) && nrow(FE[[i]]) == 0 && nrow(RE[[i]]) == 0) {
                stop("No 'para': ", pp, " found in 'model': ", mod, call. = FALSE)
            }
            if(nrow(FE[[i]]) >= 1)  FE[[i]]$model <- names(x$summary)[i]
            if(nrow(RE[[i]]) >= 1)  RE[[i]]$model <- names(x$summary)[i]
        }

        FE <-  do.call(rbind, FE)
        RE <-  do.call(rbind, RE)

        x$summary <- list("summary" =
                              list("RE" = RE,
                                   "FE" = FE,
                                   "tot_n" = NA,
                                   "convergence" = NA,
                                   "Satt_NA" = NA,
                                   "CI_NA" = NA)
                          )
    }

    if("model_selected" %in% names(res)) {
        x$model_selected <- res$model_selected
        x$model_direction <- res$model_direction
        x$LRT_alpha <- res$LRT_alpha
    }
    class(x) <- append("plcp_sim_summary", class(x))
    x
}

#' @rdname summary.plcp_sim
#' @param model Indicates which model that should be returned.
#' Default is \code{NULL} which return results from all model formulas. Can also be a \code{character} matching the
#' names used in \code{\link{sim_formula_compare}}.
#' @param model_selection indicates if the summary should be based on a LRT model selection strategy. Default is \code{NULL},
#' which returns all models, if \code{FW} or \code{BW} a forward or backward model selection strategy is used, see \emph{Details}.
#' @param LRT_alpha Indicates the alpha level used if doing LRT model comparisons.
#' @method summary plcp_sim_formula_compare
#' @export
summary.plcp_sim_formula_compare <- function(object, model = NULL, alpha = 0.05, model_selection = NULL, LRT_alpha = 0.1, para = NULL, ...) {

    if(is.null(model_selection)) {
        if(is.null(model)) {
            summary.plcp_sim(object,
                             alpha = alpha,
                             para = para,
                             ...)
        } else {
            if(is.numeric(model)) model <- names(object$res)[model]
            if(!model %in% names(object$res)) stop("No 'model' named: ", model, call. = FALSE)
            object$res <- object$res[model]
            summary.plcp_sim(object,
                             alpha = alpha,
                             para = para,
                             ...)
        }
    } else if(model_selection %in% c("FW", "BW")) {
        x <- do_model_selection(object,
                                direction = model_selection,
                                alpha = LRT_alpha)
        summary.plcp_sim(x,
                         alpha = alpha,
                         para = para,
                         ...)

    }

}


summarize_RE <- function(res, theta) {

    d <- res$RE
    parms <- unique(d$parameter)
    tmp <- vector("list", length(parms))
    for(i in seq_along(parms)) {
        para <- parms[[i]]
        vcov <- d[d$parameter == para, "vcov"]
        if(para %in% names(theta)) {
            theta_i <- theta[[para]]
        } else theta_i <- NA
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
summarize_FE <- function(res, theta, alpha, df_bw = NULL) {
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
        if(is.null(df_bw)) {
            tmp_df_bw <- unique(d[ind, "df_bw"])
        } else if(is.list(df_bw)) {
            tmp_df_bw <- df_bw[[para]]
        }

        if(any(!is.na(pval))) {
            Satt_NA <- mean(is.na(df))
        } else Satt_NA <- NA

       # para <- i
        #theta_i <- theta[[i]]
        pvals_bw <- get_p_val_df(t = estimate/se,
                                 df = tmp_df_bw,
                                 parameter = para,
                                 test = unique(d[!is.na(d$df_bw), "parameter"]))
        if(para %in% names(theta)) {
            theta_i <- theta[[para]]
        } else theta_i <- NA
        tmp[[i]] <- data.frame(
                    parameter = para,
                    M_est = mean(estimate),
                    theta = theta_i,
                    M_se = mean(se),
                    SD_est = sd(estimate),
                    Power = mean(get_cover(estimate, se, alpha = alpha)),
                    Power_bw = mean(pvals_bw < alpha),
                    Power_satt = mean(pval < alpha, na.rm=TRUE),
                    Satt_NA = Satt_NA
                )


    }
    FE <- do.call(rbind, tmp)

    FE
}
summarize_CI <- function(res, theta = NULL) {
    ## TODO: fix so theta matches fit$tests
    d <- res$FE
    parms <- unique(res$FE$parameter)

    tmp <- vector("list", length(parms))
    #thetas <- theta
    for(i in seq_along(parms)) {

        para <- parms[[i]]
        ind <- d$parameter == para
        theta_i <- theta[[para]]
        if(is.null(theta_i)) {
            theta_i <- mean(d[ind, "estimate"])
        }
        CI_lwr <- d[ind, "CI_lwr"]
        CI_upr <- d[ind, "CI_upr"]
        CI_wald_lwr <- d[ind, "CI_wald_lwr"]
        CI_wald_upr <- d[ind, "CI_wald_upr"]

        tmp[[i]] <- data.frame(
                    parameter = para,
                    CI_cover = mean(CI_lwr < theta_i &
                                         CI_upr > theta_i, na.rm = TRUE),
                    CI_Wald_cover = mean(CI_wald_lwr < theta_i &
                                             CI_wald_upr > theta_i, na.rm=TRUE),
                    CI_NA = mean(is.na(CI_lwr)))
    }
    CI_cov <- do.call(rbind, tmp)

    CI_cov
}

## Extract random effect thetas

#' @export
get_RE_thetas <- function(paras, ...) {
    UseMethod("get_RE_thetas")
}
get_RE_thetas.plcp_nested <- function(paras, ...) {
    list(
            "subject_intercept" = paras$sigma_subject_intercept^2,
            "subject_slope" = paras$sigma_subject_slope^2,
            "cluster_intercept" = paras$sigma_cluster_intercept^2,
            "cluster_slope" = paras$sigma_cluster_slope^2,
            "error" = paras$sigma_error^2,
            "cor_subject" = paras$cor_subject,
            "cor_cluster" = paras$cor_cluster
    )
}
get_RE_thetas.plcp_crossed <- function(paras, ...) {
    list(
        "subject_intercept" = paras$sigma_subject_intercept^2,
            "subject_slope" = paras$sigma_subject_slope^2,
            "cluster_intercept" = paras$sigma_cluster_intercept^2,
            "cluster_intercept_tx" = paras$sigma_cluster_intercept_tx^2,
            "cluster_slope" = paras$sigma_cluster_slope^2,
            "cluster_slope_tx" = paras$sigma_cluster_slope_tx^2,
            "error" = paras$sigma_error^2,
            "cor_subject" = paras$cor_subject,
            "cor_cluster_intercept_slope" = paras$cor_cluster_intercept_slope,
            "cor_cluster_intercept_intercept_tx" = paras$cor_cluster_intercept_intercept_tx,
            "cor_cluster_intercept_slope_tx" = paras$cor_cluster_intercept_slope_tx,
            "cor_cluster_intercept_tx_slope_tx" = paras$cor_cluster_intercept_tx_slope_tx,
            "cor_cluster_slope_intercept_tx" = paras$cor_cluster_slope_intercept_tx,
            "cor_cluster_slope_slope_tx" = paras$cor_cluster_slope_slope_tx
        )

}

## extract fixed effect thetas
#' @export
get_FE_thetas <- function(paras, ...) {
    UseMethod("get_FE_thetas")
}
get_FE_thetas.default <- function(paras, ...) {
    list(
        "(Intercept)" = paras$fixed_intercept,
        "treatment" = 0,
        "time" = paras$fixed_slope,
        "time:treatment" = get_slope_diff(paras) / paras$T_end
    )
}
summarize_convergence <- function(paras, convergence) {
    UseMethod("summarize_convergence")
}
summarize_convergence.default <- function(paras, convergence) {
    mean(convergence)
}
summary_.plcp_sim  <- function(res, paras, alpha, df_bw = NULL, ...) {
    RE_params <- get_RE_thetas(paras)

    RE <- summarize_RE(res, theta = RE_params)

    # falseconv
    false_conv <-
        res$RE[res$RE$parameter == "cluster_slope", "vcov"]
    false_conv <- mean(is_approx(false_conv, 0))

    theta <- get_FE_thetas(paras, ...)

    # support both variants of time * treatment
    t_b_t <- res$FE$parameter
    check_tbt <- t_b_t %in% c("time:treatment", "treatment:time")
    if(any(check_tbt)) {
        t_b_t <- unique(t_b_t[check_tbt])
        names(theta)[4] <- t_b_t
    }

    FE <- summarize_FE(res = res,
                       theta = theta,
                       alpha = alpha,
                       df_bw = df_bw)

    if(all(is.na(res$FE$pval))) {
        Satt_NA <- 0
    } else {
        ind <- which(!is.na(FE$Power_satt))[1]
        Satt_NA <- FE[ind, "Satt_NA"]
    }

    FE$Satt_NA <- NULL

    ## TODO: update to make compatible with fit$test
    ##       update so function is agnostic to if time:treatment or treatment:time
    CI_NA <- 0

    if ("CI_lwr" %in% colnames(res$FE)) {
        CI_cov <- summarize_CI(res, theta)
        FE$CI_Cover <- CI_cov$CI_cover
        i <- which.min(CI_cov$CI_Wald_cover >= 0)
        CI_NA <- CI_cov[i, "CI_NA"]
        FE$CI_Wald_cover <- CI_cov$CI_Wald_cover
    }

    # FE_eff <- c("(Intercept)",
    #             "treatment",
    #             "time",
    #             "time:treatment")
    #
    # FE_eff <- FE_eff[FE_eff %in% FE$parameter]
    # ind <- vapply(FE_eff, function(i)
    #                     which(FE$parameter == i), numeric(1))
    # FE <- FE[ind,]

    convergence <- summarize_convergence(paras, res$convergence)

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


.check_para_in_mod <- function(eff_names, para) {
    x <- lapply(seq_along(eff_names), function(i) {
        mod <- names(eff_names[i])
        any(eff_names[[mod]] %in% para[[mod]])
    } )

    unlist(x)
}

#' Summarize simulations based on a combination of multiple parameter values
#'
#' @param object A multiple simulation object created with
#' \code{\link{simulate.plcp_multi}}
#' @param para The name of the fixed or random effect that should be summarized.
#' @param model Specifies which model that should be summarized. Accepts either
#' a \code{character} with the name used in \code{\link{sim_formula_compare}}, or
#' an \code{integer} value.
#' @param alpha Indicates the significance level. Default is 0.05 (two-tailed),
#' one-tailed tests are not yet implemented.
#' @param model_selection Indicates if model selection should be performed. If \code{NULL} (default),
#' all models are returned, if \code{FW} or \code{BW} model selection is performed using LRT, and the result
#' is based on the selected model from each simulation. See \code{\link{summary.plcp_sim}} for more information.
#' @param LRT_alpha Indicates the alpha level used when comparing models during model selection.
#' @param ... Optional arguments.
#'
#' @method summary plcp_multi_sim
#'
#' @return A \code{list} with class \code{plcp_multi_sim_summary}. It can be coursed to a \code{data.frame},
#' using \code{as.data.frame}. Each row summarizes one of the parameter combinations used in the simulation.
#' In addition to the setup parameter values, it contains the following columns:
#' \itemize{
#'  \item \code{parameter} is the name of the coefficient
#'  \item \code{M_est} is the mean of the estimates taken over all the simulations.
#'  \item \code{theta} is the population parameter values specified with \code{study_parameters}
#'  \item \code{M_se} is the mean estimated standard error taken over all the simulations.
#'  \item \code{SD_est} is the empirical standard error; i.e. the standard
#'  deviation of the distribution of the generated estimates
#'  \item \code{power} is the empirical power of the Wald Z test, i.e. the proportion
#'  of simulated p-values < alpha
#'  \item \code{power_satt} is the empirical power of the Wald \emph{t} test using
#'   Satterthwaite's degree of freedom approximation.
#'  \item \code{satt_NA} is the proportion of Satterthwaite's approximations that failed.
#'  \item \code{prop_zero} is the proportion of the simulated estimates that
#'  are zero; only shown for random effects.
#' }
#'
#' @export
#'
summary.plcp_multi_sim <- function(object,
                                   para = "time:treatment",
                                   model = NULL,
                                   alpha = 0.05,
                                   model_selection = NULL,
                                   LRT_alpha = 0.1,
                                   ...) {
    res <- object
    mod_names <- names(res[[1]][[1]])
    mod_n <- length(mod_names)
    if (is.null(model_selection)) {
        if(!is.null(model) && is.character(model) && !model %in% mod_names) {
            stop("Incorrect 'model', no model named: ", model, call. = FALSE)
        }

        if(is.numeric(model)) {
            if(model > length(mod_names)) stop("Numeric argument 'model' is too large.", call. = FALSE)
            model <- mod_names[model]
            mod_n <- 1
        }
        if(!is.null(model) && is.character(model)) {
            if(is.list(para)) {
                if(!model %in% names(para)) stop("'model' not found in 'para'", call. = FALSE)
                para <- para[[model]]
            }

            if(is.numeric(para)) stop("'para' can't be a numeric value", call. = FALSE)
            mod_n <- 1
        }

        if(is.null(model) && is.list(para)) {
            if(!all(mod_names %in% names(para))) stop("At least one of the model names in 'para' does not exist.", call. = FALSE)
        }

    }



    # model selection
    if(!is.null(model_selection)) {
        for(i in seq_along(res)) {
            res[[i]] <- do_model_selection(res[[i]],
                                           direction = model_selection,
                                           alpha = LRT_alpha)
            model <- "model_selection"
        }
    }

    # check 'paras'
    RE_names <- lapply(res[[1]]$res, function(x) unique(x$RE$parameter))
    FE_names <- lapply(res[[1]]$res, function(x) unique(x$FE$parameter))

    if(!is.null(model)) {
        RE_names <- RE_names[[model]]
        FE_names <- FE_names[[model]]
    }

    if(length(para) == 1 && is.character(para)) {
        if(para %in% unique(unlist(RE_names))) {
            type <- "random"
        } else if(para %in% unique(unlist(FE_names))) {
            type <- "fixed"
        } else stop("No 'para' named: ", para, call. = FALSE)
    } else if(is.list(para)) {
        if(length(para) != mod_n) stop("When 'para' is a list it ",
                                       "must be the same length as the number of models: ",
                                       mod_n)
        if(any(.check_para_in_mod(RE_names, para))) {
            type <- "random"
        } else if(any(.check_para_in_mod(FE_names, para))) {
            type <- "fixed"
        } else {
            stop("At least one of the 'para'(s) was not found", call. = FALSE)
        }
    }

    # summarize
    if (type == "fixed") {
        out <- lapply(res, function(x) summary(x,
                                               para = para,
                                               model = model)$summary[[1]]$FE)
    } else if (type == "random") {
        out <- lapply(res, function(x) summary(x,
                                               para = para,
                                               model = model)$summary[[1]]$RE)
    }

    out <- as.data.frame(do.call(rbind, out))

    paras <- attr(object, "paras")

    paras <- prepare_multi_setup(paras)$out_dense
    paras <- as.data.frame(paras)
    paras <- paras[rep(1:nrow(paras), each = mod_n), ]
    out <- cbind(paras, as.data.frame(out))
    out <- out[order(out$model), ]
    class(out) <- append("plcp_multi_sim_summary", class(out))
    attr(out, "type") <- type
    attr(out, "model") <- model
    attr(out, "used_models") <- unique(out$model)
    attr(out, "used_tests") <- unique(out$parameter)
    attr(out, "nsim") <- res[[1]]$nsim
    out
}

#' Convert a multi-sim summary object to a tidy data.frame
#'
#' @param x Object with class \code{plcp_multi_sim_summary}.
#' @param ... Not used
#'
#' @return a \code{data.frame} with one row for each simulation.
#' Columns include the simulation study parameters and the results.
#' @export
as.data.frame.plcp_multi_sim_summary <- function(x, ...) {
    attr(x, "dense")
}

#' @export
as.data.frame.plcp_sim_summary <- function(x, ...) {

    RE <- lapply(seq_along(x$summary), function(i) {
        d <- x$summary[[i]]
        tmp <- d$RE
        tmp$model <- names(x$summary)[i]
        tmp$type <- "random"
        tmp
    })

}

#' Print method for \code{summary.plcp_multi_sim}-objects
#' @param x An object of class \code{plcp_multi_sim_summary}
#' @param add_cols \code{character} vector; indicates the names of the
#' additional columns that should be added to the output. Intended use case is
#' when you want to add some of the setup parameters, this print method
#' is not smart enough to figure out which parameters you are investigating.
#' @param bias \code{logical}; indicates if parameter bias should be printed.
#' @param power \code{logical}; indicates if empirical power should be printed.
#' @param estimates \code{logical}; indicates if the parameter estimates should be printed.
#' @param digits number of significant digits.
#' @param ... Optional arguments.
#' @method print plcp_multi_sim_summary
#' @export
print.plcp_multi_sim_summary <- function(x,
                                         add_cols = NULL,
                                         bias = TRUE,
                                         power = TRUE,
                                         estimates = TRUE,
                                         digits = 2, ...) {

     if(!is.null(add_cols) & !all(add_cols %in% colnames(x))) {
         not_found <- which(!add_cols %in% colnames(x))
         stop("No column called: '",
              paste(add_cols[not_found], sep = "','"),
              "'",
              call. = FALSE)
    }
    model <- attr(x, "model")
    type <- attr(x, "type")
    nsim <- attr(x, "nsim")
    out <- as.data.frame.data.frame(x)
    if(estimates) est_cols <- c("M_est", "theta", "M_se", "SD_est") else est_cols <- NULL
    if(power) power_cols <- c("Power" ,"Power_bw" ,"Power_satt") else power_cols <- NULL
    print_cols <- c("model", add_cols, est_cols, power_cols, c("CI_Cover", "CI_Wald_cover"))
    out <- out[, colnames(x) %in% print_cols]
    para <- paste(as.character(unique(attr(x, "used_tests"))), collapse = "', '")
    out <- out[, !colnames(out) %in% "parameter", ]
    mod_col <- colnames(out) %in% "model"
    if(length(attr(x, "used_models")) > 1) {
        out <- cbind(out[, mod_col, drop = FALSE],
                   out[, !mod_col])
    } else out <- out[, !mod_col]
    out[t(do.call(rbind, lapply(out, is.nan)))] <- "."
    out[t(do.call(rbind, lapply(out, is.na)))] <- "."

    if(is.null(model)) model <- "All"
    cat("Model: '", model, "' | Type: '", type, "' | Parameter(s): '", para, "'\n", sep = "")
    cat("---\n")
    print(out, digits = digits, row.names = FALSE)
    cat("---\n")
    cat("nsim: ", nsim, "|", sum(!colnames(x) %in% print_cols), "columns not shown\n")
    if(any(x$is_NA > 0)) warning("Some simulations had NA estimates that was removed.", call. = FALSE)
    invisible(out)
}

## random effect
summary_fixed.plcp_multi_sim <- function(res, para, model, alpha) {

    out <- res$res[[model]]$FE

    # for convenience catch reversed interaction
    mod_paras <- unique(out$parameter)
    TbT <- c("time:treatment", "treatment:time")
    ind <- TbT %in% mod_paras
    if(any(ind) & para %in% TbT) {
        TbT_model <- TbT[TbT %in% mod_paras]
        para <- TbT_model
    }

    theta <- switch(
        para,
        "(Intercept)" = res$paras$fixed_intercept,
        "treatment" = 0,
        "time" = res$paras$fixed_slope,
        "time:treatment" = get_slope_diff(res$paras) / res$paras$T_end,
        "treatment:time" = get_slope_diff(res$paras) / res$paras$T_end,
        NA
    )

    if(all(out$parameter != para)) stop("Argument 'para': ", "'", para, "' is not a valid parameter.", call. = FALSE)
    out <- out[out$parameter == para,]


    est <- mean(out$estimate)
    se_est <- mean(out$se)
    se_hat <- sd(out$estimate)

    power <- mean(get_cover(out$estimate,
                            out$se,
                            alpha = alpha))
    power_bw <- get_p_val_df(t = out$estimate/out$se,
                             df = out$df_bw,
                             parameter = para,
                             test = para)
    Satt_NA <-  mean(is.na(out$df))

    ret <- with(
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
            Power_satt = mean(out$pval < alpha, na.rm = FALSE),
            Satt_NA = Satt_NA

        )


    )
    if ("CI_lwr" %in% colnames(out)) {
        CI_cov <- summarize_CI(res$res[[model]], theta)
        CI_cov <- CI_cov[CI_cov$parameter == para, ]
        ret$CI_Cover <- CI_cov$CI_cover
        ret$CI_Wald_cover <- CI_cov$CI_Wald_cover
        ret$CI_NA <- CI_cov$CI_NA
    }

    ret
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

get_p_val_df <- function(t, df, parameter, test = "time:treatmen") {
    if (unique(parameter) %in% test) {
        p <- (1 - pt(abs(t), df)) * 2
    } else
        p <- NA

    p
}




