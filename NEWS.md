

# Changes in version 0.5.0

## Minor changes
* `simulate(..., seed = 1234)` no longer ignores the seed argument. If cores > 1 the "L'Ecuyer-CMRG" RNG is used via `parallel::clusterSetRNGStream`.

## Bug fixes
* `Summary(model = ...)` now correctly handles a vector of model names, and will give more informative error messages.

# Changes in version 0.4.0

This version substantially improves the `simulate` method.

## New features
* The simulate method is now much more flexible. New features include:
   - Compare more than 2 model formulas (#2).
   - Apply a data transformation during simulation; `sim_formula(..., data_transform = func)`.
     As an example `transform_to_posttest` is included.
   - Choose which parameters to test; `sim_formula(..., test = "treatment")`
   - Fit OLS models. If a model formula is supplied that contain no random effects the
   model is fit using OLS `lm()`. If this is combined with the `transform_to_posttest` the
   longitudinal model can be compared to a cross-sectional model, e.g. ANCOVA. 
   - Investigate LRT model selection of the random effects. Nested random effect models
   can be tested using LRT, and results from the "best" model is returned. The log-likelihood 
   is saved during each simulation, so the model selection can be done as a
   post-processing step; `summary.plcp_sim(..., model_selection = "FW", LRT_alpha = 0.25)`.
   
## Breaking canges
* `simulate(formula = x)` must now be created using the new functions `sim_formula`, or
`sim_formula_compare`, and can no longer be a named list or a character vector. 

## Bug fixes
* `summary.plcp_sim()` now show fixed effect `theta`s in the correct order, thanks to
GitHub user Johnzav888 (#10).

# Changes in version 0.3.0

## New features
* More flexible effect size specification. This version adds support for:
   - unstandardized effect sizes, e.g. `effect_size = 5`,
   - and Cohen's d effect sizes that are standardized 
   using either the pre- or posttest SD, or the random slope SD,
   e.g. `effect_size = cohend(0.5, "posttest_sd")`

## Other changes
* Support for lmerTest 3.0.

# Changes in version 0.2.0
## New features
* Analytical power calculations now support using Satterthwaite's degrees of
  freedom approximation.
* `Simulate.plcp` will now automatically create lme4 formulas if none is
  supplied, see `?create_lmer_formula`.
* You can now choose what alpha level to use.
* Treat cluster sizes as a random variable, `uneqal_clusters` now accepts
  a function indicating the distribution of cluster sizes, via the new argument
  `func`, e.g. `rpois` or `rnorm` could be used to draw cluster sizes.
* Expected power for designs with parameters that are random variables,
  can be calculated by averaging over multiple realizations, using the
  argument `R`.
* Support for parallel computations on Microsoft Windows, and in GUIs/interactive
  environments, using `parallel::makeCluster` (PSOCK). Forking is still used for
  non-interactive Unix environments.
  
## Improvements
* Calculations of the variance of the treatment effect is now much faster for
  designs with unequal clusters and/or missing data, when cluster sizes are
  large. The calculations now use the much faster implementation used by lme4.
* Cleaner print-methods for `plcp_multi`-objects.
* Multiple power calculations can now be performed in parallel, via the
  argument `cores`.
* `simulate.plcp_multi` now have more options for saving intermediate results.
* `print.plcp_multi_power` now has better support for subsetting via either [], 
head(), or subset().

## Breaking changes
* `icc_pre_subject` is now defined as `(u_0^2 + v_0^2) / (u_0^2 + v_0^2 + error^2)`,
instead of `(u_0^2) / (u_0^2 + v_0^2 + error^2)`. This would be the subject-level ICC, 
if there's no random slopes, i.e. correlation between time points for the same subject.
* `study_parameters()`: 0 and NA now means different things. If 0 is passed, the parameters
is kept in the model, if you want to remove it specify it as NA instead.
* `study_parameters()`: is now less flexible, but more robust. Previously a large
combination if raw and relative parameters could be combined, and the individual 
parameters was solved for. To make the function less bug prone and easier to maintain,
it is now only possible to specify the cluster-level variance components as relative values,
if the other parameters as passed as raw inputs. 

## Bug fixes and minor changes
* Output from `simulate_data()` now includes a column `y_c` that contains the full outcome vector,
without missing values added. This makes it easy to compare the complete and incomplete
data set, e.g. via `simulate()`.
* `simulate()` new argument `batch_progress` enables showing progress when doing
  multiple simulations.
* Bugfix for `summary.plcp_sim` where the wrong % convergence was calculated.
* Simulation function now accepts lme4 formulas containing "||".
* The cluster-level intercept variance is now also set to zero in the control
  group, when a partially nested design is requested.
* Fix incorrect error message from `study_parameters` when
  `icc_cluster_pre = NULL` and all inputs are standardized.
* Fix bug that would cause all slopes to be zero when `var_ratio` argument was
  passed a vector of values including a 0, e.g. `var_ratio = c(0, 0.1, 0.2)`.
* Bugfix for multi-sim objects that caused the wrong class the be used for,
  e.g. `res[[1]]$paras`, and thus the single simulation would not print
  correctly.
* Results from multi-sim objects can now be summarized for all random effects
  in the model.
* More support for summarizing random effects from partially nested formulas,
  e.g. `cluster_intercept` and `cluster_slope` is now correctly extracted from
  `(0 + treatment + treatment:time || cluster)`.
* When Satterthwaite's method fails the between clusters/subjects *dfs*
  are used to calculate *p*-values.
* `Power.plcp_multi` is now exported.
* `get_power.plcp_multi` now shows a progress bar.
* Fix bug that caused dropout to be wrong when one condition had 0 dropout, and
`deterministic_dropout = FALSE`.

# Changes in version 0.1.0
First release.
