# powerlmm 0.2.0.9000


## New features
* Analytical power calculations now support using Satterthwaite's degrees of freedom approximation. 
* Simulate.plcp will now automatically create lmer formulas if none is supplied,
see ?create_lmer_formula
* You can now choose what alpha level to use. 
* Threat cluster sizes as a random variable. `uneqal_clusters` now accepts
 a function indicating the distribution of cluster sizes, via the new argument
`func`, e.g. `rpois` or `rnorm` could be used to draw cluster sizes.
* Expected power for designs with parameters that are random variables, 
can be calculated by averaging over multiple realizations, using the argument `R`. 
* Support for parallel computions on Microsoft Windows, and in GUIs/interactive environments,
using `parallel::makeCluster`. Forking is still used for non-interactive Unix environments.
  
## Improvements
* Calculations of the variance of treatment effects is now much faster for designs with 
unequal clusters and/or missing data, when cluster sizes are large. The calculations now
use the much faster implementation used by lme4.
* Cleaner print-methods for `plcp_multi`-objects

## Minor changes and fixes
* Fixed some typos in the documentation
* New arguments to simulate.plcp_multi; save_folder, etc
* 'simulate.pclp' new argument 'batch_progress' enables showing progress when doing
multiple simulations.
* Fixed bug in summary.plcp_sim, wrong % convergence calculated. 
* Simulation function now accepts lmer-formulas containing "||"
* Bug fix, the cluster-level intercept variance is now also set to zero in control group,
for partially nested designs.
* Fixed incorrect error message from study_parameters when 'icc_cluster_pre = NULL' an all inputs are
standardized
* Fixed a bug that would cause all slopes to be zero when 'var_ratio' argument was 
passed a vector of values including a 0, e.g. var_ratio = c(0, 0.1, 0.2)
* Fixed a bug for multi_sim-objects that caused the wrong class the be used for, e.g. 'res[[1]]$paras', and thus
used the wrong print method.
* Results from mult-sim objects can now be summarized for all random effects in the model.
* More support for summarizing random effects from partially nested formulas,
e.g. 'cluster_intercept' and 'cluster_slope' is now correctly extracted from
 '(0 + treatment + treatment:time || cluster)' 
* Simulation, if Satterthwaite's method fails DF are calculated as total n2 - 2.
* Power.plcp_multi is now exported
* Get_power.plcp_multi now shows a progress bar
