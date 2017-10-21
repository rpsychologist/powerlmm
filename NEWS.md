# powerlmm 0.2.0.9000


## New features
* Simulate.plcp will now automatically create lmer formulas if none is supplied,
see ?create_lmer_formula
* Analytical power calculations now support using Satterthwaite's degrees of freedom approximation.

## Improvements
* Calculations of the variance of treatment effects is now much faster for designs with 
unequal clusters and/or missing data, when cluster sizes are large. The calculations now
use the much faster implementation used by lme4.

## Minor changes and fixes
* Fixed some typos in the documentation
* New arguments to simulate.plcp_multi; save_folder, etc
* 'simulate.pclp' new argument 'batch_progress' enables showing progress when doing
multiple simulations.
* Fixed bug in summary.plcp_sim, wrong % convergence calculated. 
* Simulation function now accepts lmer-formulas containing "||"
* Bug fix, the cluster-level intercept variance is now also set to zero in control group,
for partially nested designs.
* Fixed incorrect error message from study_parameters when icc_cluster_pre = NULL an all inputs are
standardized
* Fixed a bug that would cause all slopes to be zero when 'var_ratio' argument was 
passed a vector of values including a 0, e.g. var_ratio = c(0, 0.1, 0.2)
