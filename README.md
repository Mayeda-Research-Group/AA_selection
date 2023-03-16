# AA_selection

This repo contains R scripts for the paper named Estimating dementia incidence in older Asian Americans and Pacific Islanders in California: an application of inverse odds of selection weights. 

1. `imputation_in_rpgeh.R` performs imputation for some demographic and health-related variables in the RPGEH dataset. 

2. Scripts for data preparation: 

* `2.1_chis_combine_harmonize.R` combines CHIS data from 2005, 2007, and 2009, and harmonizes them with RPGEH.
* `2.2_rpgeh_imputed_harmonize.R` harmonizes RPGEH data with CHIS.
* `2.3_table1.R` prepares tables of summary statistics pre- and post-imputation. 

3. Covariate balance and weight development: 

* `3.1_covariate_balance_and_weight_generation.R` checks covariate balance between unweighted RPGEH and CHIS, and then generates weights for RPGEH and check covariate balance iteratively. Many of the results are exploratory, and are NOT used for estimation of incidence rates.
* `3.2_summary_pre_post_weighting.R` generates covariate balance plots by ethnicity before and after applying weight and summary statistics for propensity scores and weights. 

4. Scripts to calculate IR with weights and bootstrap for confidence interval: 

* `4.1_bootstrap_IR_prep.R` prepares the dataset for IR calculation. 
* `4.2_bootstrap_IR.R` sets up the IR calculation to be run on the Hoffman2 cluster. 
* `4.3_submission_script.sh` is the job submission script on the cluster. 
* `4.4_bootstrap_processing.R` cleans the output. 

5. `sensitivity_age_in_weighting_models.R` contains sensitivity analysis on how to include age in the GLM weighting model. 

Auxiliary scripts include: 

* `functions_summary.R` contains functions called in `table1.R`. 
* `function_IR_calc.R` contains helper functions for IR calculation. 

