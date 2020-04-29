
print_per_treatment <- function(n, width = 0, n2 = FALSE, hanging = 19) {

    x <- lapply(seq_along(n), function(i) {
        print_per_treatment_(i, x = n, n2 = n2)$lab
    })
    x <- format(x, width = width)
    x <- paste(x, " (", names(n), ")", sep ="")
    collapse <- paste0("\n", paste(rep(" ", hanging), collapse = ""))
    x <- paste(unlist(x), collapse = collapse)
    x
}
print_per_treatment_ <- function(i, x, n2 = FALSE) {
    name <- names(x)[i]
    x <- x[[i]]
    if(n2 & length(unique(x)) == 1) {
        if(attr(x, "func")) {
            x_num <- NA
            x_lab <- paste(unique(x))
        } else {
            x_num <- unique(x)
            if(length(x) > 1) {
                x_lab <- paste(unique(x),"x", length(x))
            } else {
                x_lab <- paste(unique(x))
            }

        }
    } else if(length(unique(x)) == 1) {
        x_lab <- x
        x_num <- x
    } else {
        x_lab <- x
        x_num <- NA
    }
    list("lab" = paste(paste(unlist(x_lab), collapse = ", "), sep =""),
         "num" = x_num)
}

deparse_n2 <- function(n2) {
    n2_attr <- attr(n2, "func")
    if(!is.null(n2_attr) && (n2_attr != "manual")) {
        n2 <- deparse(n2_attr)
        attr(n2, "func") <- TRUE
    } else attr(n2, "func") <- FALSE

    n2

}
truncate_n2 <- function(n2) {
    n2 <- gsub("^.*::", "", n2)
    n2 <- strtrim(n2, 20)

}
prepare_print_n2 <- function(x) {
    n2 <- get_n2(x)
    n2$treatment <- deparse_n2(n2$treatment)
    n2_attr <- attr(n2$control, "func")
    if(attr(n2$control, "per_treatment")) {
        n2$control <- deparse_n2(n2$control)

    } else if(!is.null(n2_attr) && (n2_attr != "manual")) {
        n2$control <- "-"
        attr(n2$control, "func") <- FALSE
    } else {
        n2$control <- deparse_n2(n2$control)
    }


    n2
}
prepare_print_plcp <- function(x, two_level = FALSE, crossed = FALSE, hanging = 19) {
    n1 <- x$n1
    n2 <- prepare_print_n2(x)
    n3 <- get_n3(x)
    if(crossed) n3 <- n3["total"]/2
    tot_n <- get_tot_n(x)
    width <- max(nchar(print_per_treatment_(1, n2, n2 = TRUE)$lab),
                 nchar(print_per_treatment_(2, n2, n2 = TRUE)$lab),
                 nchar(print_per_treatment_(3, tot_n)$lab))
    if(two_level) width <- max(vapply(tot_n, nchar, numeric(1)))
    n2 <- print_per_treatment(n2,
                              width = width,
                              n2 = TRUE,
                              hanging = hanging)
    n3 <- print_per_treatment(n3,
                              width = width,
                              hanging = hanging)

    tot_n <- print_per_treatment(tot_n,
                                 width = width,
                                 hanging = hanging)

    icc_slope <- round(get_ICC_slope(x), 2)
    var_ratio <- round(get_var_ratio(x), 2)
    icc_pre_clusters <- round(get_ICC_pre_clusters(x), 2)
    icc_pre_subjects <- round(get_ICC_pre_subjects(x), 2)

    effect <- get_effect_size(x)
    effect_label <- ifelse(effect$standardizer == "raw", "raw", "Cohen's d")
    if(x$partially_nested) {
        effect_label <- ifelse(effect_label == "Cohen's d",
                               paste(effect_label, " [SD: ", effect$standardizer,
                                     ", ", effect$treatment, "]", sep = ""),
                               effect_label)
    } else {
        effect_label <- ifelse(effect_label == "Cohen's d",
                               paste(effect_label, " [SD: ", effect$standardizer, "]", sep = ""),
                               effect_label)
    }

    effect <- paste(effect$ES, " (", effect_label,")", sep = "")

    gd <- get_dropout(x)
    gd$time <- format(gd$time, nsmall = 0, digits = 3, width = 2)
    gd$control <-  format(gd$control*100, nsmall = 0, digits = 0, width = 2)
    gd$treatment <- format(gd$treatment*100, nsmall = 0, digits = 0, width = 2)
    colnames(gd) <- c("time", "%, control", "%, treatment")
    gd <- print_per_treatment(gd,
                              hanging = hanging)
                              
    method <- ifelse(crossed,
        "Study setup (three-level, crossed)",
        "Study setup (three-level, nested)")

      res <- structure(list(n1 = n1,
            n2 = n2,
            n3 = n3,
            total_n = tot_n,
            dropout = gd,
            icc_pre_subjects = icc_pre_subjects,
            icc_pre_clusters = icc_pre_clusters,
            icc_slope = icc_slope,
            var_ratio = var_ratio,
            `effect_size` = effect,
            method = method),
        class = "power.htest")
    if (crossed) {
        icc_slope_crossed <- round(
            get_ICC_slope(x, interaction = TRUE),
            2)
        icc_pre_clusters_crossed <- round(
            get_ICC_pre_clusters(x, interaction = TRUE),
            2)
        res$icc_slope <- paste0(icc_slope, " (crossed: ", icc_slope_crossed, ")")
        res$icc_pre_clusters <- paste0(icc_pre_clusters, " (crossed: ", icc_pre_clusters_crossed, ")")

    } else {
  
    }

    attr(res, "width") <- width
    res
}
prepare_print_plcp_2lvl <- function(x) {
    res <- prepare_print_plcp(x, two_level = TRUE)
    if(!is.list(x$dropout)) res$dropout <- "No missing data"
    res$method <- "Study setup (two-level)"
    res$icc_slope <- NULL
    res$icc_pre_clusters <- NULL
    res$n2 <- res$total_n
    res$n3 <- NULL
    res$total_n <- NULL

    res
}
prepare_print_plcp_hurdle_2lvl <- function(x) {
    res <- prepare_print_plcp(x, two_level = TRUE, hanging = 18)
    if(!is.list(x$dropout)) res$dropout <- "No missing data"

    if(x$marginal) marginal <- " marginal" else marginal <- NULL
    res$method <- paste0("Study setup (two-level,", marginal, " hurdle-", x$family, ")")
    res$icc_slope <- NULL
    res$icc_pre_clusters <- NULL
    res$n2 <- res$total_n
    res$n3 <- NULL
    res$total_n <- NULL
    res$icc_pre_subjects <- NULL
    res$var_ratio <- NULL


    ES_cont <- paste0("Continuous (multiplicative): ", x$RR_cont)
    res$effect_size <- paste0("Zeros (OR): ", x$OR_hu,
                              "\n", paste0(rep(" ", 18), collapse = ""),
                              ES_cont)

    names(res)[names(res) == "effect_size_2"] <- ""

    # attrs <- attributes(res)
    #
    # # reorder
    # res <- structure(res[c("n1", "n2", "dropout", "family", "effect_size", "method")],
    #                  class = "power.htest")
    # attr(res, "width") <- attrs$width
    #



    res

}

prepare_print_plcp_3lvl <- function(x, crossed = FALSE) {
    res <- prepare_print_plcp(x, crossed = crossed)
    if(!is.list(x$dropout)) res$dropout <- "No missing data"
    if(x$partially_nested) res$method <- "Study setup (three-level, partially nested)"
    if(is.unequal_clusters(x$n2)) {
        if(is.per_treatment(x$n2)) {
            if(class(x$n2[[1]]$treatment[[1]]) == "plcp_unequal_clusters") {
                n2 <- x$n2[[1]]$treatment[[1]]()
                n2_attr <- attr(n2, "func")
            } else {
                n2 <- x$n2[[1]]$control[[1]]()
                n2_attr <- attr(n2, "func")
            }
        } else {
            n2 <- x$n2[[1]]()
            n2_attr <- attr(n2, "func")
        }

        if(!is.null(n2_attr) &
           n2_attr != "manual") res$note <- "n2 is randomly sampled"
        if(!is.null(n2_attr) &
           n2_attr != "manual" &
           !is.per_treatment(x$n2)) attr(res, "same_dist") <- TRUE
    }
    res
}

#' Print method for three-level \code{study_parameters}-objects
#' @param x An object of class \code{plcp_3lvl}.
#' @param ... Optional arguments.
#' @method print plcp_3lvl
#' @export
print.plcp_3lvl <- function(x, ...) {
    res <- prepare_print_plcp_3lvl(x)
    print(res, digits = 2)
    if(!is.null(attr(res, "same_dist"))) message("The same random 'n2' sample is used for both treatment groups,\n'per_treatment' is needed to sample both groups independently.")
}

print.plcp_crossed <- function(x, ...) {
    res <- prepare_print_plcp_3lvl(x, crossed = TRUE)
    print(res, digits = 2)
}

#' Print method for two-level \code{study_parameters}-objects
#' @param x An object of class \code{plcp_2lvl}.
#' @param ... Optional arguments.
#' @method print plcp_2lvl
#' @export
print.plcp_2lvl <- function(x, ...) {
   res <- prepare_print_plcp_2lvl(x)
   print(res, digits = 2, ...)
}

#' Print method for two-part/hurdle \code{study_parameters}-objects
#' @param x An object of class \code{plcp_hurdle}.
#' @param ... Optional arguments.
#' @method print plcp_2lvl
#' @export
print.plcp_hurdle <- function(x, ...) {
    res <- prepare_print_plcp_hurdle_2lvl(x)
    print(res, digits = 2, ...)
}