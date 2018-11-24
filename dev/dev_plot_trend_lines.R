
p <- study_parameters(n1 = 11,
                      n2 = 5,
                      n3 = 1,
                      icc_pre_subject = 0.5,
                      var_ratio = 0.03,
                      icc_slope = 0.05,
                      cor_subject = -0.8,
                      cor_cluster = -0.5,
                      effect_size = cohend(-0.5))

d <- simulate_data(p)

ggplot(d, aes(time, mu2, group = subject)) + geom_line(alpha =.3) +
    stat_summary(fun.y = "mean", geom = "line", aes(group = treatment, color = treatment), size = 2)


ggplot(d, aes(time, mu3, group = subject, color = treatment)) + geom_line(alpha =.3) +
    stat_summary(fun.y = "mean", geom = "line", aes(group = treatment), size = 2)


ggplot(d,aes(time, mu2, group = subject)) +
    geom_point(aes(y=y)) +
    geom_line(alpha =.3) +
    #stat_summary(fun.y = "mean", geom = "line", aes(group = treatment, color = treatment), size = 2) +
    facet_wrap(~subject)
