# powerlmm 0.1.1.9000


## New features
* Simulate.plcp will now automatically create lmer formulas if none is supplied,
see ?create_lmer_formula

## Minor changes
* Fixed some types in the documentation
* New arguments to simulate.plcp_multi; save_folder, etc
* 'simulate.pclp' new argument 'batch_progress' enables showing progress when doing
multiple simulations.
* Fixed bug in summary.plcp_sim, wrong % convergance calculated. 
* Simulation function now accepts lmer-formulas containing "||"

