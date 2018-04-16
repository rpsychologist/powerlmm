<!-- README.md is generated from README.Rmd. Please edit that file -->
powerlmm
========

[![Travis-CI Build Status](https://travis-ci.org/rpsychologist/powerlmm.svg?branch=master)](https://travis-ci.org/rpsychologist/powerlmm) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/powerlmm)](https://cran.r-project.org/package=powerlmm)

Power Analysis for Longitudinal Multilevel/Linear Mixed-Effects Models.

Overview
--------

The purpose of `powerlmm` is to help design longitudinal treatment studies (parallel groups), with or without higher-level clustering (e.g. longitudinally clustered by therapists, groups, or physician), and missing data. The main features of the package are:

-   Longitudinal two- and three-level (nested) linear mixed-effects models, and partially nested designs.
-   Random slopes at the subject- and cluster-level.
-   Missing data.
-   Unbalanced designs (both unequal cluster sizes, and treatment groups).
-   Design effect, and estimated type I error when the third-level is ignored.
-   Fast analytical power calculations for all designs.
-   Power for small samples sizes using Satterthwaite's degrees of freedom approximation.
-   Explore bias, Type I errors and model misspecification using convenient simulation methods.

Installation
------------

`powerlmm` can be installed from CRAN and GitHub.

``` r
# CRAN, version 0.2.0
install.packages("powerlmm")

# GitHub, dev version
devtools::install_github("rpsychologist/powerlmm", build_vignettes = TRUE)
```

Example usage
-------------

This is an example of setting up a three-level longitudinal model with random slopes at both the subject- and cluster-level, with different missing data patterns per treatment arm. Relative standardized inputs are used, but unstandardized raw parameters values can also be used.

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
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

p
#> 
#>      Study setup (three-level) 
#> 
#>               n1 = 11
#>               n2 = 10 x 5 (treatment)
#>                    10 x 5 (control)
#>               n3 = 5      (treatment)
#>                    5      (control)
#>                    10     (total)
#>          total_n = 50     (treatment)
#>                    50     (control)
#>                    100    (total)
#>          dropout =  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 (time)
#>                     0,  0,  1,  3,  6,  9, 12, 16, 20, 25, 30 (%, control)
#>                     0,  0,  1,  2,  4,  5,  8, 10, 13, 17, 20 (%, treatment)
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>      effect_size = -0.8 (Cohen's d [SD: pretest_SD])
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
#>               n2 = 10 x 5 (treatment)
#>                    10 x 5 (control)
#>               n3 = 5      (treatment)
#>                    5      (control)
#>                    10     (total)
#>          total_n = 50     (control)
#>                    50     (treatment)
#>                    100    (total)
#>          dropout =  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 (time)
#>                     0,  0,  1,  3,  6,  9, 12, 16, 20, 25, 30 (%, control)
#>                     0,  0,  1,  2,  4,  5,  8, 10, 13, 17, 20 (%, treatment)
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>      effect_size = -0.8 (Cohen's d [SD: pretest_SD])
#>               df = 7.947873
#>            alpha = 0.05
#>            power = 68%
```

### Unequal cluster sizes

Unequal cluster sizes is also supported, the cluster sizes can either be random (sampled), or the marginal distribution can be specified.

``` r
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(2, 3, 5, 20),
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

get_power(p)
#> 
#>      Power Analyis for Longitudinal Linear Mixed-Effects Models (three-level)
#>                   with missing data and unbalanced designs 
#> 
#>               n1 = 11
#>               n2 = 2, 3, 5, 20 (treatment)
#>                    2, 3, 5, 20 (control)
#>               n3 = 4           (treatment)
#>                    4           (control)
#>                    8           (total)
#>          total_n = 30          (control)
#>                    30          (treatment)
#>                    60          (total)
#>          dropout = No missing data
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>      effect_size = -0.8 (Cohen's d [SD: pretest_SD])
#>               df = 6
#>            alpha = 0.05
#>            power = 44%
```

Cluster sizes follow a Poisson distribution

``` r
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(func = rpois(5, 5)), # sample from Poisson
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))

get_power(p, R = 100, progress = FALSE) # expected power by averaging over R realizations
#> 
#>      Power Analyis for Longitudinal Linear Mixed-Effects Models (three-level)
#>                   with missing data and unbalanced designs 
#> 
#>               n1 = 11
#>               n2 = rpois(5, 5) (treatment)
#>                    -           (control)
#>               n3 = 5           (treatment)
#>                    5           (control)
#>                    10          (total)
#>          total_n = 24.89       (control)
#>                    24.89       (treatment)
#>                    49.78       (total)
#>          dropout = No missing data
#> icc_pre_subjects = 0.5
#> icc_pre_clusters = 0
#>        icc_slope = 0.05
#>        var_ratio = 0.02
#>      effect_size = -0.8 (Cohen's d [SD: pretest_SD])
#>               df = 8
#>            alpha = 0.05
#>            power = 48% (MCSE: 1%)
#> 
#> NOTE: n2 is randomly sampled. Values are the mean from R = 100 realizations.
```

### Convenience functions

Several convenience functions are also included, e.g. for creating power curves.

``` r
x <- get_power_table(p, 
                     n2 = 5:10, 
                     n3 = c(4, 8, 12), 
                     effect_size = cohend(c(0.5, 0.8), standardizer = "pretest_SD"))
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
