#' @importFrom matrixStats rowMeans2 rowMedians
#' @export
marginalize <- function(object, R, ...) {
    UseMethod("marginalize")
}
