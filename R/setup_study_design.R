# Study design ---------------------------------------------------

#' @export
study_design <- function(nested = TRUE,
                         levels = 3,
                         time_form = "linear",
                         custom = FALSE) {

    if (time_form != "linear") {
        message("Argument 'time_form' is currently ignored.", call. = FALSE)
        time_form <- "linear"
    }

    # check inputs
    stopifnot(is.logical(nested))
    stopifnot(levels %in% 2:3)
    stopifnot(time_form == "linear")

    args <- list(
        nested = nested,
        levels = levels,
        time_form = time_form
        )

    if (custom) {
        class(args) <- "plcp_design_custom"
    } else if (nested) {
        class(args) <- "plcp_design_nested"

    } else {
        class(args) <- "plcp_design_crossed"
    }
    class(args) <- append(class(args), "plcp_design")

    args
}


print.plcp_design <- function(x, ...) {
    if (x$levels == 3) {
        levels <- paste0("3 ", ifelse(x$nested, "(nested)", "(crossed)"))
    } else if (x$levels == 2) {
        levels <- 2
    }
    res <- structure(list(levels = levels,
        time_form = x$time_form,
        method = "Study design"),
    class = "power.htest")

    print(res)
}