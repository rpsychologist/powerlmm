## New formula syntax

list(list("label" = "y ~ time + (1 | cluster)",
     data_transform = function(),
     package = "lme4"
     )
     )

# or
list("label1" = list(formula  ="y ~ time + (1 | cluster)",
                    data_transform = function(),
                    package = "lme4"
                    ),
     "label2" = list(formula  ="y ~ pre + treatment + (1 | cluster)",
                     data_transform = function(),
                     package = "lme4"
                     )
     )

## backward compat
## must work with
f <- list("correct" = "y ~ treatment * time",
          "wrong" = "y ~ treatment * time + (1 + time | subject)")

## better to use helper functions; transparent arguments
list("label1" = sim_formula(formula  ="y ~ time + (1 | cluster)",
                     data_transform = function(),
                     package = "lme4"
),
"label2" = sim_formula(formula  ="y ~ pre + treatment + (1 | cluster)",
                data_transform = function(),
                package = "lme4"
                )
)

# advantage of also using helper instead of outer list()
# custom print method for formula object

# need seperate func to combine sim_formula
# compare_sim_formulas(sim_formula, sim_formula)
# can also add LRT argument to compare_sim_formula()
#
# compare_sim_formula(sim_formula(...), sim_formula(...), LRT = TRUE)
#
# also add alpha arg?


# some stepwise sim function? conditional removal of cluster effects?
# perhaps best to view it as a post processing function, if LRT result are
# available?


# TODO --------------------------------------------------------------------
# * enable LRT comparison of elements of compare_sim_formula
