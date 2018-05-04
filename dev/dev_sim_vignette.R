# Same model as we used in the previous example
p <- study_parameters(n1 = 11,
                      n2 = 40,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      var_ratio = 0.02,
                      dropout = c(0, dropout_weibull(0.3, 1)),
                      effect_size = cohend(0))

# simulation formulas
# analyze as a posttest only 2-level model
f_pt <- sim_formula("y ~ treatment",
                    test = "treatment",
                    data_transform = transform_to_posttest)

f_pt_pre <- sim_formula("y ~ treatment + pretest",
                    test = "treatment",
                    data_transform = transform_to_posttest)

# analyze as 3-level longitudinal
f_lt <- sim_formula("y ~ time*treatment +
                         (1 + time | subject)")
# constrained
f_lt2 <- sim_formula("y ~ time + time:treatment +
                         (1 + time | subject)")

f <- sim_formula_compare("posttest" = f_pt,
                         "post_covariate" = f_pt_pre,
                         "longitudinal" = f_lt,
                         "longi_constrained" = f_lt2)

res <- simulate(p,
                formula = f,
                nsim = 15,
                cores = 15,
                satterthwaite = TRUE)

# extract different para per model
# 'para' is not needed, the default is to print all model parameters.
summary(res, model = NULL,
        para = list("posttest" = "treatment",
                    "post_covariate" = "treatment",
                    "longitudinal" = "time:treatment",
                    "longi_constrained" = "time:treatment"))


summary(res, model = NULL,
        para = list("posttest" = "treatment",
                    "post_covariate" = "treatment",
                    "longitudinal" = "time:treatment",
                    "longi_constrained" = "time:treatment"))

summary(res, model = 1, para = "treatment")

summary(res, para = "subject_slope")
summary(res, para = "subject_slope", model = 3)

# should give error
summary(res, model = "posttest",
        para = list("posttest" = "treatment",
                    "post_covariate" = "treatment",
                    "longitudinal" = "time:treatment",
                    "longi_constrained" = "time:treatment"))

res1 <- simulate(p,
                formula = f_pt,
                nsim = 5,
                cores = 1,
                satterthwaite = TRUE)




