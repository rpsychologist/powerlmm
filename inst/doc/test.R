## ------------------------------------------------------------------------
library(powerlmm)
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
    eval = NOT_CRAN
    )

NOT_CRAN

## ------------------------------------------------------------------------
Sys.getenv("NOT_CRAN")

## ------------------------------------------------------------------------
if(NOT_CRAN) TRUE else stop("print IS_CRAN")

