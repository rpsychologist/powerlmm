study_parameters.plcp_design_crossed <- function(design,
                                     n1,
                                     n2,
                                     n3 = 1,
                                     T_end = NULL,
                                     fixed_intercept = 0L,
                                     fixed_tx = 0L,
                                     fixed_slope = 0L,
                                     sigma_subject_intercept = NULL,
                                     sigma_subject_slope = NULL,
                                     sigma_cluster_intercept = NULL,
                                     sigma_cluster_slope = NULL,
                                     sigma_cluster_intercept_tx = NULL,
                                     sigma_cluster_slope_tx = NULL,
                                     sigma_error = 10,
                                     cor_subject = 0L,
                                     cor_cluster_intercept_slope = 0L,
                                     cor_cluster_intercept_intercept_tx = 0L,
                                     cor_cluster_intercept_slope_tx = 0L,
                                     cor_cluster_slope_intercept_tx = 0L,
                                     cor_cluster_slope_slope_tx = 0L,
                                     cor_cluster_intercept_tx_slope_tx = 0L,
                                     cor_within = 0L,
                                     var_ratio = NULL,
                                     icc_slope = NULL,
                                     icc_pre_subject = NULL,
                                     icc_pre_cluster = NULL,
                                     effect_size = 0L,
                                     cohend = NULL,
                                     dropout = 0L,
                                     deterministic_dropout = TRUE,
                                     ...) {


    # deprecated Cohen's d
    if(!is.null(cohend)) {
        effect_size <- cohend(cohend, standardizer = "pretest_SD", treatment = "control")
    }

    if(is.per_treatment(n3)) stop("the number of clusters can't differ per preatment in a crossed design.", call. = FALSE)


    args <- list(
        n1 = n1,
        n2 = n2,
        n3 = n3,
        T_end = T_end,
        fixed_intercept = fixed_intercept,
        fixed_tx = fixed_tx,
        fixed_slope = fixed_slope,
        sigma_subject_intercept = sigma_subject_intercept,
        sigma_subject_slope = sigma_subject_slope,
        sigma_cluster_intercept = sigma_cluster_intercept,
        sigma_cluster_slope = sigma_cluster_slope,
        sigma_cluster_intercept_tx = sigma_cluster_intercept_tx,
        sigma_cluster_slope_tx = sigma_cluster_slope_tx,
        sigma_error = sigma_error,
        cor_subject = cor_subject,
        cor_cluster_intercept_slope = cor_cluster_intercept_slope,
        cor_cluster_intercept_intercept_tx = cor_cluster_intercept_intercept_tx,
        cor_cluster_intercept_slope_tx = cor_cluster_intercept_slope_tx,
        cor_cluster_slope_intercept_tx = cor_cluster_slope_intercept_tx,
        cor_cluster_slope_slope_tx = cor_cluster_slope_slope_tx,
        cor_cluster_intercept_tx_slope_tx = cor_cluster_intercept_tx_slope_tx,
        cor_within = cor_within,
        var_ratio = var_ratio,
        icc_slope = icc_slope,
        icc_slope_tx = icc_slope_tx,
        icc_pre_subject = icc_pre_subject,
        icc_pre_cluster = icc_pre_cluster,
        icc_pre_cluster_tx = icc_pre_cluster_tx,
        effect_size = effect_size,
        dropout = dropout,
        deterministic_dropout = deterministic_dropout
    )
    save_call <- args

    tmp_args <- args[!vapply(args, is.null, logical(1))]



    # ## Default NA
    # if(is.null(icc_pre_subject) & is.null(sigma_subject_intercept)) {
    #     tmp_args$icc_pre_subject <- NA
    # }
    # if(is.null(var_ratio) & (is.null(sigma_subject_slope) || is.na(sigma_subject_slope)) & is.null(sigma_cluster_slope)) {
    #     tmp_args$var_ratio <- NA
    # }
    # if(is.null(icc_pre_cluster) & is.null(sigma_cluster_intercept)) {
    #     tmp_args$icc_pre_cluster <- NA
    # }
    # if(is.null(icc_slope) & is.null(sigma_cluster_slope)) {
    #     tmp_args$icc_slope <- NA
    # }
     tmp <- expand.grid(tmp_args)
    # ## --
    #
    # ## icc_cluster > icc_subject
    # if(!is.null(icc_pre_cluster) & !is.null(icc_pre_subject)) {
    #     if(any((tmp$icc_pre_cluster > tmp$icc_pre_subject), na.rm=TRUE)) {
    #         stop("'icc_pre_cluster' can't be larger than 'icc_pre_subject'", call. = FALSE)
    #     }
    # }
    #
    # # check cluster slope variance exists when var ratio is NA or zero.
    # if(!is.null(tmp$var_ratio)) {
    #     if(is.na(tmp$var_ratio) && any(tmp$icc_slope >= 0, na.rm=TRUE)) {
    #         stop("'icc_slope' can't be >= 0 when 'var_ratio' or 'sigma_subject_slope' is NA", call. = FALSE)
    #     }
    #     if((any(tmp$var_ratio == 0, na.rm=TRUE) | any(sigma_subject_slope == 0)) &&
    #        any(tmp$icc_slope > 0, na.rm=TRUE)) {
    #         stop("'icc_slope' can't be > 0 when 'var_ratio' or 'sigma_subject_slope' is 0", call. = FALSE)
    #     }
    # }

    ## Default T_end
    if(is.null(args$T_end)) tmp$T_end <- tmp$n1 - 1


    ## Solve raw values

    # # Solve subject_intercept
    # if(is.null(sigma_subject_intercept)) {
    #
    #     icc_cluster <- tmp$icc_pre_cluster
    #     icc_cluster[is.na(icc_cluster)] <- 0
    #
    #     tmp$sigma_subject_intercept <-
    #         sqrt(
    #             tmp$sigma_error^2 /
    #                 (1-(tmp$icc_pre_subject)) *
    #                 (tmp$icc_pre_subject - icc_cluster)
    #         )
    # }
    #
    # # solve subject_slope
    # if(is.null(sigma_subject_slope)) {
    #
    #     icc <- tmp$icc_slope
    #     icc[is.na(icc) | is.null(icc)] <- 0
    #     tmp$sigma_subject_slope <- sqrt(tmp$var_ratio *
    #                                         tmp$sigma_error^2 * (1-icc))
    #
    # }
    # # Solve cluster_intercept
    # if(is.null(sigma_cluster_intercept)) {
    #     # from icc_pre_subject
    #     if(is.null(sigma_subject_intercept)) {
    #         tmp$sigma_cluster_intercept <-
    #             sqrt(
    #                 tmp$sigma_error^2 /
    #                     (1-(tmp$icc_pre_subject)) *
    #                     tmp$icc_pre_cluster
    #             )
    #         # from sigma_subject_intercept
    #     } else {
    #         v0 <- with(tmp, (icc_pre_cluster*(sigma_error^2 + sigma_subject_intercept^2))
    #                    /(1-icc_pre_cluster))
    #         tmp$sigma_cluster_intercept <- sqrt(v0)
    #     }
    #
    # }
    #
    # # Solve cluster_slope
    # if(is.null(sigma_cluster_slope)) {
    #
    #     # check if NA
    #     if(is.null(icc_slope) || all(is.na(icc_slope)) ) {
    #
    #         tmp$sigma_cluster_slope <- NA
    #     } else {
    #         # solve from icc_slope and var_ratio
    #         if(is.null(sigma_subject_slope)) {
    #             v1 <- with(tmp, var_ratio * sigma_error^2 * icc_slope)
    #             tmp$sigma_cluster_slope <- sqrt(v1)
    #             # solve from subject_slope and icc_slope
    #         } else {
    #             x <- with(tmp, sigma_subject_slope^2/(1-icc_slope))
    #             v1 <- x - tmp$sigma_subject_slope^2
    #             v1 <- sqrt(v1)
    #             tmp$sigma_cluster_slope <- v1
    #
    #         }
    #     }
    # }
    #

    # keep cols
    cols <- which(colnames(tmp) %in% c("icc_slope",
                                       "var_ratio",
                                       "icc_pre_cluster",
                                       "icc_pre_subject"))
    cols <- colnames(tmp)[cols]
    paras <- tmp[, !(names(tmp) %in% cols)]


    # Single or multi?
    if((is.data.frame(paras) & nrow(paras) == 1)) {
        paras <- as.list(paras)
    }

    paras$design <- "plcp_design_crossed"
    # avoid problems with prepare_paras()
    paras$partially_nested <- FALSE
    paras$custom_model <- FALSE
    # if(is.data.frame(paras)) {
    #     class(paras) <- append(c("plcp_multi"), class(paras))
    # } else class(paras) <- append(c("plcp"), class(paras))
    #
    # # Default cor_*
    # if(is.null(paras$cor_cluster)) paras$cor_cluster <- cor_cluster
    # if(is.null(paras$cor_subject)) paras$cor_subject <- cor_subject

    # # Classes
    # if(all(is.na(paras$sigma_cluster_slope)) &
    #    all(is.na(paras$sigma_cluster_intercept))) {
    #     class(paras) <- append(class(paras), c("plcp_2lvl"))
    # } else if (all(!is.na(paras$sigma_cluster_slope)) |
    #            all(!is.na(paras$sigma_cluster_intercept))) {
    #     class(paras) <- append(class(paras), c("plcp_3lvl"))
    # } else {
    #     class(paras) <- append(class(paras), c("plcp_mixed"))
    # }

    class(paras) <- append(class(paras), c("plcp","plcp_crossed"))

    attr(paras, "call") <- save_call
    paras

}

