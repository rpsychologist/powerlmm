# Study design ---------------------------------------------------

#' @export
study_design <- function(nested = TRUE,
                         family = "gaussian",
                         levels = 3,
                         groups = 2,
                         time_form = "linear",
                         custom = FALSE) {

    if(groups != 2) {
        message("Argument 'groups' is currently ignored.", call. = FALSE)
        groups <- 2
    }
    if(time_form != "linear") {
        message("Argument 'time_form' is currently ignored.", call. = FALSE)
        time_form <- "linear"
    }
    if(any(family %in% c("hurdle", "two-part"))
       & !nested) stop("'crossed' designs are not yet implemented for hurdle models", call. = FALSE)
    if(any(family %in% c("hurdle", "two-part"))
       & levels != 2) {
        warning("3 level designs are not yet implemented for hurdle models", call. = FALSE)
        levels <- 2
    }



    families <- c("gaussian",
                  "binomial", "poisson",
                  "gamma", "lognormal",
                  "hurdle", "two-part")
    if(!family %in% families) stop("Not a supported 'family'", call. = FALSE)

    # check inputs
    stopifnot(is.logical(nested))
    stopifnot(levels %in% 2:3)
    stopifnot(time_form == "linear")


     if(family %in% c("two-part", "hurdle")) {
         family_class <- "hurdle"
     } else {
         family_class <- NULL
     }

    args <- list(nested = nested,
                 levels = levels,
                 groups = groups,
                 time_form = time_form,
                 family = family)

    if (custom) {
        class(args) <- "plcp_design_custom"
    } else if(nested) {
        class(args) <- paste(c("plcp_design", family_class, "nested"),
                             collapse = "_")

    } else {
        class(args) <- paste(c("plcp_design", family_class, "crossed"),
                             collapse = "_")
    }
    class(args) <- append(class(args), "plcp_design")

    args
}


print.plcp_design <- function(x, ...) {
    if(x$levels == 3) {
        levels <- paste0("3 ", ifelse(x$nested, "(nested)", "(crossed)"))
    } else if(x$levels == 2) {
        levels <- 2
    }
    res <- structure(list(levels = levels,
                          groups = x$groups,
                          family = x$family,
                          time_form = x$time_form,
                          method = "Study design"),
                     class = "power.htest")

    print(res)
}
