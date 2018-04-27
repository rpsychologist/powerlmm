## single
models <- list(list("label" = "m0",
                    "ll" = -10,
                    "df" = 1),
               list("label" = "m1",
                    "ll" = -8,
                    "df" = 2),
               list("label" = "m2",
                    "ll" = -8,
                    "df" = 3),
               list("label" = "m3",
                    "ll"= -6,
                    "df" = 4)
)
step_bw.plcp_sim(models) # should pick 'm3'
step_fw.plcp_sim(models) # should pick 'm1'


## multiple
models <- list(list("label" = "m0",
                    "ll" = c(-10, -10),
                    "df" = c(1,1)),
               list("label" = "m1",
                    "ll" = c(-10, -8),
                    "df" = c(2,2)),
               list("label" = "m2",
                    "ll" = c(-10, -10),
                    "df" = c(3,3)),
               list("label" = "m3",
                    "ll"= c(-10, -10),
                    "df" = c(4,4))
)
step_bw.plcp_sim(models) # m0, m1
step_fw.plcp_sim(models) # m0, m1

#
models <- list(list("label" = "m0",
                    "ll" = c(-10, -10),
                    "df" = c(1,1)),
               list("label" = "m1",
                    "ll" = c(-10, -8),
                    "df" = c(2,2)),
               list("label" = "m2",
                    "ll" = c(-6, -10),
                    "df" = c(3,3)),
               list("label" = "m3",
                    "ll"= c(-10, -10),
                    "df" = c(4,4))
)
step_bw.plcp_sim(models) # m2, m1
step_fw.plcp_sim(models) # m0, m1



models <- list(list("label" = "m0",
                    "ll" = c(-10, -10),
                    "df" = c(1, 1)),
               list("label" = "m1",
                    "ll" = c(-8, -8),
                    "df" = c(2, 2)),
               list("label" = "m2",
                    "ll" = c(-6, -6),
                    "df" = c(3, 3)),
               list("label" = "m3",
                    "ll"= c(-4, -5),
                    "df" = c(4, 4))
)
step_bw.plcp_sim(models) # m3, m2
step_fw.plcp_sim(models) # m3, m2





## sim
p <- study_parameters(n1 = 11,
                      n2 = unequal_clusters(5,4,20,20,30),
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      icc_slope = 0.05,
                      var_ratio = 0.02,
                      dropout = dropout_weibull(0.3, 0.3),
                      effect_size = cohend(c(0, 0.8)))

f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
f1 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
f2 <- sim_formula("y ~ time * treatment + (1 + time | subject) + (0 + time | cluster)")
f <- compare_sim_formulas("m0" = f0, "m1" = f1, "m2" = f2)


res <- simulate(p, formula = f, nsim = 100, satterthwaite = TRUE, cores = 10, CI = FALSE)
summary(res)

x <-  summary(res)


summary(res, model_selection = "FW")
summary(res, model = "m2")
