


# Binomial ----------------------------------------------------------------
p_bin <- study_parameters(design = study_design(family = "binomial"),
                          n1 = 11,
                          n2 = 25,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.02,
                          effect_size = log(1),
                          sigma_error = 1)

m_bin <- marginalize(p_bin)


# Poisson ----------------------------------------------------------------
p_pois <- study_parameters(design = study_design(family = "poisson"),
                          n1 = 11,
                          n2 = 25,
                          icc_pre_subject = 0.5,
                          var_ratio = 0.02,
                          effect_size = log(1),
                          sigma_error = 1)

m_pois <- marginalize(p_pois)

# lognormal ----------------------------------------------------------------
p_ln <- study_parameters(design = study_design(family = "lognormal"),
                           n1 = 11,
                           n2 = 25,
                           icc_pre_subject = 0.5,
                           var_ratio = 0.02,
                           effect_size = log(1),
                           sigma_error = 1)

m_ln <- marginalize(p_ln)

# Gamma ----------------------------------------------------------------
p_gamma <- study_parameters(design = study_design(family = "gamma"),
                         n1 = 11,
                         n2 = 25,
                         icc_pre_subject = 0.5,
                         var_ratio = 0.02,
                         effect_size = log(1),
                         shape = 1.5,
                         sigma_error = 1)

m_gamma <- marginalize(p_gamma)
