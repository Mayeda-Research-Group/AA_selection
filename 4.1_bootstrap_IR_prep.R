# set up data and various objects in preparation of bootstrapping
# including the concatenated chis and rpgeh data (with 40 imps), 
# glm formulas by ethnicity for weight generation, 
# age categories and calculated pys data for IR calculation. 

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "twang", "boot", "openxlsx"
       # "mice", rlang",
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "openxlsx", "lmtest", "mitools", "ggpubr", "patchwork", "RColorBrewer"
)

options(scipen = 999, digits = 8)

# set up paths and data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_RPGEH <- "Asian_Americans_dementia_data/analysis_data_tables/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"

rpgeh_tte <- read_sas(paste0(path_to_box, path_to_RPGEH, 
                             "aa_adrd_cardiometabolic_tte.sas7bdat"))
rpgeh_harm_imp <- readRDS(paste0(path_to_box, path_to_datasets, 
                                 "rpgeh_imp_harmonized.RDS")) 
chis_harm <- readRDS(paste0(path_to_box, path_to_datasets, 
                            "chis_2005_to_2009_harmonized.RDS"))

# stack data for GLM and weight generation ----
# load the function that labels and refactors all the variables of interest
source(paste0(path_to_box, 
              "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
              "cleaned_scripts/functions_summary.R"))

chis_rpgeh <- chis_harm %>% 
  mutate(imp = 0) %>% 
  select(-INTVLANG) %>% 
  add_row(rpgeh_harm_imp) %>% 
  t1_relabel()

# RPGEH tte data for py calculation ----
# load the function that calculates person-years and case contributions
source(paste0(path_to_box, 
              "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
              "cleaned_scripts/function_IR_calc.R"))

rpgeh_dem <- rpgeh_tte %>% 
  select( 
    # ID and tte variables
    SUBJID, SURVEY_AGE, 
    MAIN_DEM_V1_END_AGE, #MAIN_DEM_V1_END_DEM_AGE, 
    MAIN_DEM_V1_END_DEM_FLAG) %>% 
  # reformat ID variable into string
  mutate(ID = as.character(SUBJID), .before = "SURVEY_AGE") %>% 
  select(-SUBJID) # drop old ID variable

# calculate person-years and case contribution by age category
rpgeh_dem_pys <- calculate_pys(
  data = rpgeh_dem, 
  age_cat = c(59, 65, 70, 75, 80, 85, 120), 
  fu_start = "SURVEY_AGE", fu_end = "MAIN_DEM_V1_END_AGE", 
  event_flag = "MAIN_DEM_V1_END_DEM_FLAG")

# output is: 
pys_data <- rpgeh_dem_pys$data # pys and cases by age category
age_cat_labels <- rpgeh_dem_pys$age_cat_labels # lavels for the age categories
# pys and case contribution are the same over the bootstrapping process

# GLM formulas used for each ethnicity for weight generation ----
fmls <- c(
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_hhsize_3 + H_health_3 + H_hyp", 
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_hhsize_3 + H_hyp + H_bmi + H_retired + H_health_3", 
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_hhsize_3 + H_hyp + H_health_3",
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_health_3 + H_hyp",
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_health_3 + H_hyp",
  "in_rpgeh ~ H_age + H_female + H_edu_4 + H_retired + H_health_3 + H_hyp",
  "in_rpgeh ~ H_age + H_female + H_edu_4 + H_work + H_health_3 + H_bmi + H_hyp",
  "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_smk + H_health_3"
)

names(fmls) <- c("Chinese", "Filipino", "Japanese", "Korean", "Pacific Islander",
                 "South Asian", "Vietnamese", "Multiple Ethnicities")


# save the various objects for running on Hoffman ----
# save(chis_rpgeh, pys_data, age_cat_labels, fmls,
#      file = paste0(path_to_box, path_to_datasets, "bootstrap_prep_data.RData")
#      )

# preliminary IR calculation ----
# we can calculate unweighted age-specific pys and case contributions 
# for each ethnicity since they do not depend on imputation or the glm models.
# this helps us identify age groups within each ethnicity with fewer than 
# 5 events, which will be suppressed in the manuscript results

pys_data_w_ethn <- chis_rpgeh %>% 
  filter(imp == 1) %>% 
  select(ID, H_ethn) %>% 
  left_join(pys_data, by = "ID") 

unw_ncases_ethn <- pys_data_w_ethn %>% 
  group_by(H_ethn) %>% 
  summarize(
    across(
      starts_with("case"),
      .fns = sum,
      .names = "n_{.col}"
    )
  )

# write.xlsx(unw_ncases_ethn,
#            paste0(path_to_box, 
#                   "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
#                   "cleaned_scripts/output/bootstrap_results_strat/", 
#                   "unweighted_n_cases_x_ethn.xlsx"))


### keep everything above as prep before Hoffman ###
# the following is OLD code!


# bootstrapping function ----
# for one imputation, and for one bootstrap sample
# (1) weight generation 
# data: take one imputation of RPGEH and complete CHIS
# add bootstrap sampling (either write it myself, or use indices for boot())
# calculate unconditional sampling odds 
# fit final glm and calculate weights
# (2) IR calculation
# join in tte data by subject
# calculate crude and age-adjusted IRs, unweighted and with weights
# output is 4 statistics: uwt_crude, wt_crude, uwt_adj, wt_adj

#### ----

# # argument
# fml <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4 + H_hhsize_3 + H_hyp + H_health_3" 
# 
# # take Japanese as an example
# # prepare the input data for imp = 1
# # (this is prepared before the bootstrap function)
# input_data <- chis_rpgeh %>% 
#   filter(H_ethn == "Japanese", imp %in% c(0, 1))
# 
# boot_data <- input_data[indices, ]
# # boot_data <- input_data
# 
# # calculate unconditional RPGEH sampling odds
# p_rpgeh <- sum(boot_data$in_rpgeh) / sum(boot_data$smplwt)
# uncon_rpgeh_odds <- p_rpgeh / (1 - p_rpgeh)
# 
# # fit the model
# mod <- glm(as.formula(fml), data = boot_data,
#            family = binomial(link = "logit"), weights = boot_data$smplwt)
# # calculate predicted probabilities
# # boot_data$p1 <- predict.glm(mod, type = "response")
# boot_data <- boot_data %>% 
#   # generate weights
#   mutate(
#     p1 = predict.glm(mod, type = "response"), 
#     iow1 = (1 - p1) / p1,
#     sw1 = ifelse(in_rpgeh == 1, iow1 * uncon_rpgeh_odds, 1),
#     # sw1_twang = sw1 * smplwt # this is only for twang - bal.stat
#   ) %>% 
#   # we don't need CHIS data for IR calculations
#   filter(in_rpgeh == 1) %>% 
#   # join RPGEH pys data 
#   left_join(., pys_data, "ID")
# 
# calculate_ir(data = boot_data, 
#              age_cat = age_cat_labels, 
#              std_pop = c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
#              )
# 
# 
