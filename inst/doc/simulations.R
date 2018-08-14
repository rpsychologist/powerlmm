## ---- include = FALSE-------------------------------------------------------------------
library(knitr)
library(dplyr)

# set 'Sys.setenv("VIGNETTE_EVAL" = TRUE)' to print cached sims
VIGNETTE_EVAL <- identical(tolower(Sys.getenv("VIGNETTE_EVAL")), "true")

# set 'Sys.setenv("RUN_SIM" = TRUE)' re-run sims
RUN_SIM <- identical(tolower(Sys.getenv("RUN_SIM")), "true")

# avoid inline eval error
if(!VIGNETTE_EVAL) {
    round2 <- function(...) NA
} else {
    round2 <- round
}

options(width = 90)
knitr::opts_chunk$set(
    eval = VIGNETTE_EVAL,
    purl = VIGNETTE_EVAL
    )

## ---------------------------------------------------------------------------------------
#cores <- parallel::detectCores(logical = FALSE) # use all physical CPU cores
cores <- 16
nsim <- 1000

## ---- echo = FALSE----------------------------------------------------------------------
if(!RUN_SIM & file.exists("../dev/vignette_cache/sim_res.Rdata")) load("../dev/vignette_cache/sim_res.Rdata")
if(!RUN_SIM & file.exists("dev/vignette_cache/sim_res.Rdata")) load("dev/vignette_cache/sim_res.Rdata")

## ---------------------------------------------------------------------------------------
library(powerlmm)
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 6,
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = 0.05,
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                      partially_nested = TRUE,
                      effect_size = cohend(-0.8, 
                                           standardizer = "pretest_SD"))
p

## ---------------------------------------------------------------------------------------
x <- get_power(p)
get_monte_carlo_se(x, nsim = c(1000, 2000, 5000))

## ---------------------------------------------------------------------------------------
f <- sim_formula("y ~ treatment * time + 
                 (1 + time | subject) + 
                 (0 + treatment:time | cluster)")

if(RUN_SIM) {
    res <- simulate(object = p,
                    nsim = nsim,
                    formula = f,
                    satterthwaite = TRUE,
                    cores = cores,
                    save = FALSE)
    
}

summary(res)

## ---------------------------------------------------------------------------------------
p <- update(p, effect_size = 0)

f <- sim_formula_compare("correct" = "y ~ treatment * time + 
                                      (1 + time | subject) + 
                                      (0 + treatment:time | cluster)",
                         "wrong" = "y ~ treatment * time + 
                                    (1 + time | subject)")
if(RUN_SIM) {
    res2 <- simulate(object = p,
                     nsim = nsim,
                     formula = f,
                     satterthwaite = TRUE,
                     cores = cores,
                     save = FALSE)
}


summary(res2)

## ---------------------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = c(6, 12),
                      fixed_intercept = 37,
                      fixed_slope = -0.64,
                      sigma_subject_intercept = 2.8,
                      sigma_subject_slope = 0.4,
                      sigma_cluster_intercept = 0,
                      cor_subject = -0.2,
                      icc_slope = c(0.05, 0.1),
                      sigma_error = 2.6,
                      dropout = dropout_weibull(proportion = 0.3, 
                                                rate = 1/2),
                      partially_nested = TRUE,
                      effect_size = cohend(c(0.5, 0.8), 
                                           standardizer = "pretest_SD"))


f <- sim_formula("y ~ treatment * time + (1 + time | subject) + 
                 (0 + treatment:time | cluster)")
if(RUN_SIM) {
    res3 <- simulate(object = p,
                     nsim = nsim,
                     formula = f,
                     satterthwaite = TRUE,
                     cores = cores,
                     batch_progress = FALSE)
}


## ---------------------------------------------------------------------------------------
# Summarize the 'time:treatment' results 
x <- summary(res3, para = "treatment:time")

# customize what to print
print(x, add_cols = c("n3_tx", "icc_slope", "effect_size"),
      estimates = FALSE)

# Summarize the cluster-level random slope 
x <- summary(res3, para = "cluster_slope")
print(x, add_cols = c("n3_tx", "icc_slope", "effect_size"))

## ---------------------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 40,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      var_ratio = c(0, 0.01, 0.02, 0.03),
                      dropout = c(0, dropout_weibull(0.3, 1)),
                      effect_size = cohend(c(0, 0.8)))

# Formulas --------------------------------------------------------------------
# OLS "t-test"
f_PT <- sim_formula("y ~ treatment",
                    test = "treatment",
                    data_transform = transform_to_posttest)

# ANCOVA
f_PT_pre <- sim_formula("y ~ treatment + pretest",
                        test = "treatment",
                        data_transform = transform_to_posttest)

# analyze as 2-level longitudinal
f_LMM <- sim_formula("y ~ time * treatment +
                         (1 + time | subject)")

# constrain treatment differences at pretest
f_LMM_c <- sim_formula("y ~ time + time:treatment +
                         (1 + time | subject)")

# combine formulas
f <- sim_formula_compare("posttest" = f_PT,
                         "ANCOVA" = f_PT_pre,
                         "LMM" = f_LMM,
                         "LMM_c" = f_LMM_c)

# Run sim --------------------------------------------------------------------
if(RUN_SIM) {
    res4 <- simulate(p,
                     formula = f,
                     nsim = nsim,
                     cores = cores,
                     satterthwaite = TRUE,
                     batch_progress = FALSE) 
}

## ---------------------------------------------------------------------------------------
tests <-  list("posttest" = "treatment",
               "ANCOVA" = "treatment",
               "LMM" = "time:treatment",
               "LMM_c" = "time:treatment")

x <- summary(res4, para = tests)
print(head(x, 5), 
      add_cols = c("var_ratio", 
                   "effect_size"))


## ---- dpi = 150, fig.asp = 0.8, fig.width = 7, out.width = "100%", fig.align = "center", warnings = FALSE, message = FALSE----
library(ggplot2)
x$model <- factor(x$model, levels = c("posttest",
                                      "ANCOVA", 
                                      "LMM", 
                                      "LMM_c")) 

d_limits <- data.frame(effect_size = c(0), 
                       Power_satt = c(0.025, 0.075),
                       var_ratio = 0, 
                       model = 0,
                       dropout = 0)

d_hline <- data.frame(effect_size = c(0, 0.8), 
                      x = c(0.05, 0.8))

ggplot(x, aes(model, 
               Power_satt,
               group = interaction(var_ratio, dropout), 
               color = factor(var_ratio),
               linetype = factor(dropout))) + 
    geom_hline(data = d_hline,  aes(yintercept = x)) +
    geom_point() + 
    geom_blank(data = d_limits) +
    geom_line() +
    labs(y = "Power (Satterthwaite)",
         linetype = "Dropout",
         color = "Variance ratio",
         title = "Power and Type I errors",
         subtitle = 'Comparing "cross-sectional" and longitudinal models\nunder various amounts of heterogeneity and dropout') +
    facet_grid(dropout ~ effect_size, 
               scales = "free", 
               labeller = "label_both") +
    coord_flip()

## ---------------------------------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 40,
                      n3 = 3,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.03)

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f3 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (1 + time | cluster)")

f <- sim_formula_compare("subject-intercept" = f0, 
                         "subject-slope" = f1, 
                         "cluster-slope" = f2,
                         "cluster-intercept" = f3)


## ---------------------------------------------------------------------------------------
if(RUN_SIM) {
    res5 <- simulate(p, formula = f, 
                     nsim = nsim, 
                     satterthwaite = TRUE, 
                     cores = cores,
                     CI = FALSE)    
}


## ---------------------------------------------------------------------------------------
x <- summary(res5, model_selection = "FW", LRT_alpha = 0.05)
x

## ---- dpi = 150, fig.asp = 0.8, fig.width = 7, out.width = "100%", fig.align = "center", warnings = FALSE, message = FALSE----
alphas <- seq(0.01, 0.5, length.out = 50)
x <- vapply(alphas, function(a) {
    type1 <- summary(res5, model_selection = "FW", LRT_alpha = a)
    type1$summary$model_selection$FE$Power_satt[4]
    }, numeric(1))

d <- data.frame(LRT_alpha = alphas, type1 = x)
ggplot(d, aes(LRT_alpha, type1)) + 
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    labs(y = "Type I errors (time:treatment)",
         title = "Impact of the alpha level when using LRT for model selection")

## ---------------------------------------------------------------------------------------
x1 <- summary(res5, para = "time:treatment")
x1

## ---------------------------------------------------------------------------------------
# See if power is impacted
p1 <- update(p, effect_size = cohend(0.8))
if(RUN_SIM) {
    res6 <- simulate(p1, 
                     formula = f, 
                     nsim = nsim, 
                     satterthwaite = TRUE,
                     cores = cores, 
                     CI = FALSE) 
}


# we can also summary a fixed effect for convenience
x <- summary(res6, model_selection = "FW", para = "time:treatment")
print(x, verbose = FALSE)
x1 <- summary(res6, para = "time:treatment")
print(x1, verbose = FALSE)

## ---------------------------------------------------------------------------------------
# simulation formulas

# analyze as a posttest only 1-level OLS model
f0 <- sim_formula("y ~ treatment",
                 test = "treatment",
                 data_transform = transform_to_posttest)

# analyze as a posttest only 2-level model
f1 <- sim_formula("y ~ treatment + (1 | cluster)",
                 test = "treatment",
                 data_transform = transform_to_posttest)

f <- sim_formula_compare("post_ignore" = f0, 
                         "post_2lvl" = f1)

if(RUN_SIM) {
    res7 <- simulate(p,
                     formula = f,
                     nsim = nsim,
                     cores = cores,
                     satterthwaite = TRUE,
                     batch_progress = FALSE) 
}


# Type I error
print(summary(res7, para = "treatment"), verbose = FALSE)

## ---------------------------------------------------------------------------------------
# model selection
summary(res7,
        model_selection = "FW", 
        para = "treatment", 
        LRT_alpha = 0.1)

## ---- echo = FALSE----------------------------------------------------------------------
if(RUN_SIM) save(res, res2, res3, res4, res5, res6, res7, 
                 file = "../dev/vignette_cache/sim_res.Rdata")

