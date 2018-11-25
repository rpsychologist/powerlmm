

plot(m, RE = TRUE,
     type = "trend",
     outcome = "overall",
     overlay = TRUE, #
     RE_level = c(1,2), trim = c(0, 0.99),
     fixed_subject_percentiles = c(0.5, 0.5),
     fixed_subject_hu_percentiles = c(0.01)) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 100, 1000, 10000))


cowplot::plot_grid(
    plot(m, RE = TRUE, type = "trend_ridges",
     RE_level = c(1), trim = c(0, 0.999),
     fixed_subject_percentiles = c(0.5, 0.5),
     fixed_subject_hu_percentiles = c(0.5)) +
     scale_x_continuous(trans = "sqrt",
                        breaks = c(0, 1, 100, 1000, 10000), limits = c(0, 2e4)),
    plot(m, RE = TRUE, type = "trend_ridges",
         RE_level = c(1), trim = c(0, 0.99),
         fixed_subject_percentiles = c(0.95, 0.5),
         fixed_subject_hu_percentiles = c(0.5)) +
        scale_x_continuous(trans = "sqrt",
                           breaks = c(0, 1, 100, 1000, 10000), limits = c(0, 2e4))
     )



# allow multiple sd percentiles
# maybe next release
plot(m,
     RE = TRUE,
     type = "trend_ridges",
     RE_level = c(1,2), trim = c(0, 0.99),
     fixed_subject_percentiles =list(c(0.5, 0.1),
                                     c(0.5, 0.9)),
     fixed_subject_hu_percentiles = c(0.5)) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 100, 1000, 10000))
