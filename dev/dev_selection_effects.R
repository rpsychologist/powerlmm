p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      partially_nested = TRUE,
                      effect_size = 0)

trans <- function(sd) {

    function(data) {
        d <- data
        d_tx <- d[d$treatment == 1, ]
        n3 <- length(unique(d_tx$cluster))

        pre <- d_tx[d_tx$time == 0, "subject_intercept"]

        pre <- pre + rnorm(length(pre), 0, sd)
        cluster <- as.numeric(cut(pre,
                                  breaks = quantile(pre, seq(0, 1, length.out = n3 + 1)),
                                  include.lowest = TRUE,
                                  ordered = TRUE))

        d_tx$cluster <- cluster[d_tx$subject]

        d_cc <- d[d$treatment == 0, ]
        d_cc$treatment <- 0
        d_cc$cluster <- n3 + 1
        d_cc$subject <- d_cc$subject + nrow(d_tx)

        d <- rbind(d_cc, d_tx)

        d

    }

}

trans_f <- trans(sd = 20)
#f <- sim_formula("y ~ time + (1 | subject) + (1 | cluster)", data_transform = trans_f, test = "time")

f <- sim_formula("y ~ time * treatment  + (1 | subject) + (0 + treatment | cluster)", data_transform = trans_f)
f1 <- sim_formula("y ~ time * treatment  + (1 | subject)")
res <- simulate(p, nsim = 500, formula = sim_formula_compare("selection" = f, "2-lvl" = f1), cores = 16)
summary(res)

summary(res, para = "time:treatment")

# true nesting
p2 <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0.1,
                      partially_nested = TRUE,
                      effect_size = 0)

f2 <- sim_formula("y ~ time * treatment  + (1 | subject) + (0 + treatment | cluster)")
f3 <- sim_formula("y ~ time * treatment  + (1 | subject)")

res2 <- simulate(p2, nsim = 500, formula = sim_formula_compare("PN" = f2, "2lvl" = f3), cores = 16)
summary(res2)


# post-test

p2 <- study_parameters(n1 = 11,
                       n2 = 50,
                       n3 = 4,
                       icc_pre_subject = 0.5,
                       icc_pre_cluster = 0.1,
                       partially_nested = TRUE,
                       effect_size = 0)

trans_selection_post <- function(func) {
    func <- func
    function(data) {
        data <- func(data)

        transform_to_posttest(data)
    }

}
trans_f2 <- trans_selection_post(func= trans(sd = 20))


f_post1 <- sim_formula("y ~ treatment", test = "treatment", data_transform = transform_to_posttest)
f_post2 <- sim_formula("y ~ treatment + (0 + treatment | cluster)", test = "treatment", data_transform = transform_to_posttest)



f_post_selection1 <- sim_formula("y ~ treatment",
                                 test = "treatment",
                                 data_transform = trans_f2)
f_post_selection2 <- sim_formula("y ~ treatment + (0 + treatment | cluster)",
                                 test = "treatment",
                                 data_transform = trans_f2)


res3 <- simulate(p2, nsim = 500,
                 formula = sim_formula_compare("1lvl" = f_post1,
                                              "2lvl" = f_post2),
                 cores = 16)
summary(res3)

res4 <- simulate(p, nsim = 500,
                 formula = sim_formula_compare("1lvl" = f_post_selection1,
                                               "2lvl" = f_post_selection2),
                 cores = 16)
summary(res4)





# Random slope ------------------------------------------------------------
p <- study_parameters(n1 = 11,
                      n2 = 50,
                      n3 = 4,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.03,
                      cor_subject = -0.5,
                      partially_nested = TRUE,
                      effect_size = 0)


trans_f <- trans(sd = 20)
f <- sim_formula("y ~ time * treatment  + (1 + time | subject) + (0 + treatment | cluster)",
                 data_transform = trans_f)
f3 <- sim_formula("y ~ time * treatment  + (1 + time | subject) + (0 + treatment + treatment:time | cluster)",
                 data_transform = trans_f)
f1 <- sim_formula("y ~ time * treatment  + (1 + time | subject)")
res <- simulate(p,
                nsim = 500,
                formula = sim_formula_compare("selection" = f,
                                              "2-lvl" = f1,
                                              "3-lvl" = f3),
                cores = 16)
summary(res)

summary(res, para = "time:treatment")

# ES ----------------------------------------------------------------------

trans_f(simulate_data(p)) %>%
    filter(time == 10) %>%
    group_by(treatment, time) %>% summarise(mean(y)) %>%
    print(n = 20)



# plot --------------------------------------------------------------------
d <- trans(20)(simulate_data(p))
tmp <- d[d$treatment == 1 & d$time == 0, ]

plot(tmp$cluster, tmp$subject_intercept)

tmp <- d %>%
    filter(treatment == 1) %>%
    mutate(mu = subject_intercept + subject_slope * time)
tmp_agg <- tmp %>%
    group_by(cluster,time) %>%
    summarise(mu = mean(mu))

tmp %>%
    ggplot(aes(time, mu, color = cluster, group = subject)) +
    geom_line() +
    geom_line(data = tmp_agg, aes(group = cluster), color = "red", size =1)




