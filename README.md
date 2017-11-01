<!-- README.md is generated from README.Rmd. Please edit that file -->
powerlmm
========

[![Travis-CI Build Status](https://travis-ci.org/rpsychologist/powerlmm.svg?branch=master)](https://travis-ci.org/rpsychologist/powerlmm) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/powerlmm)](https://cran.r-project.org/package=powerlmm)

Power Analysis for Longitudinal Multilevel/Linear Mixed Models.

Overview
--------

The purpose of `powerlmm` is to help design longitudinal treatment studies, with or without higher-level clustering (e.g. longitudinally clustered by therapists, groups, or physician), and missing data. The main features of the package are:

-   Longitudinal two- and three-level (nested) linear mixed models, and partially nested designs.
-   Random slopes at the subject- and cluster-level.
-   Missing data.
-   Unbalanced designs (both unequal cluster sizes, and treatment groups).
-   Design effect, and estimated type I error when the third-level is ignored.
-   Fast analytical power calculations for all designs.
-   Power for small samples sizes using Satterthwaite's degrees of freedom approximation.
-   Explore bias, Type 1 errors and model misspecification using convenient simulation methods.

Installation
------------

`powerlmm` can be installed from CRAN and GitHub.

``` r
# CRAN
install.packages("powerlmm")

# GitHub
devtools::install_github("rpsychologist/powerlmm", build_vignettes = TRUE)
```

Example usage
-------------

This is an example of setting up a three-level longitudinal model with random slopes at both the subject- and cluster-level, with different missing data pattern per treatment arm. Here relative standardized inputs are used, but unstandardized raw parameters values can also be used.

``` r
library(powerlmm)
d <- per_treatment(control = dropout_weibull(0.3, 2),
              treatment = dropout_weibull(0.2, 2))
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 5,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      dropout = d,
                      cohend = -0.8)

p
#> 
#>      Study setup (three-level) 
#> 
#>               n1 = 11
#>               n2 = 10  (treatment)
#>                    10  (control)
#>               n3 = 5   (treatment)
#>                    5   (control)
#>                    10  (total)
#>          total_n = 50  (treatment)
#>                    50  (control)
#>                    100 (total)
#>          dropout =  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 (time)
#>                     0,  0,  1,  3,  6,  9, 12, 16, 20, 25, 30 (%, control)
#>                     0,  0,  1,  2,  4,  5,  8, 10, 13, 17, 20 (%, treatment)
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>           cohend = -0.8
```

``` r
plot(p)
```

![](http://rpsychologist.com/img/powerlmm/README-three-level-setup-1.png)

``` r
get_power(p, df = "satterthwaite")
#> 
#>      Power Analyis for Longitudinal Linear Mixed-Effects Models (three-level)
#>                   with missing data and unbalanced designs 
#> 
#>               n1 = 11
#>               n2 = 10  (treatment)
#>                    10  (control)
#>               n3 = 5   (treatment)
#>                    5   (control)
#>                    10  (total)
#>          total_n = 50  (treatment)
#>                    50  (control)
#>                    100 (total)
#>          dropout =  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 (time)
#>                     0,  0,  1,  3,  6,  9, 12, 16, 20, 25, 30 (%, control)
#>                     0,  0,  1,  2,  4,  5,  8, 10, 13, 17, 20 (%, treatment)
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>           cohend = -0.8
#>               df = 8.8139
#>            alpha = 0.05
#>            power = 69 %
```

Several convenience functions are also included, e.g. for creating power curves.

``` r
x <- get_power_table(p, n2 = 5:10, n3 = c(4, 8, 12), cohend = c(0.5, 0.8))
```

``` r
plot(x)
```

![](http://rpsychologist.com/img/powerlmm/README-three-level-power-curve-1.png)

Launch interactive web application
----------------------------------

The package's basic functionality is also implemented in a Shiny web application, aimed at users that are less familiar with R. Launch the application by typing

``` r
library(powerlmm)
shiny_powerlmm()
```

![](http://rpsychologist.com/img/powerlmm/README-shiny-screenshot1.png)

![](http://rpsychologist.com/img/powerlmm/README-shiny-screenshot2.png)
