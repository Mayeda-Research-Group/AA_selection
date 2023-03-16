# imputation for RPGEH
# Note: CHIS is complete

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "mice"
       # rlang"
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "openxlsx", "ggpubr", "mitools", "lmtest"
)

options(scipen = 999, digits = 8)

# set up paths and data ----
# desktop paths
# path_to_box <- "C:/Users/Staff/Box/"
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_RPGEH <- "Asian_Americans_dementia_data/analysis_data_tables/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
rpgeh_tte <- read_sas(paste0(path_to_box, path_to_RPGEH, 
                             "aa_adrd_cardiometabolic_tte.sas7bdat"))

# prepare dataset ----
pre_MI_data <- rpgeh_tte %>% 
  filter(
    # exclude dementia cases pre-survey (prevalent dementia) 
    # and cases with follow-up time = 0
    MAIN_DEM_V1_SAMPLE == 1,
    # include Asians only
    !is.na(ETHNICITY_REV), ETHNICITY_REV != 9
  )
# check
# table(pre_MI_data$ETHNICITY_REV, useNA = "ifany")

# missingness summary in RPGEH, Asians only ----
## harmonizable variables ----
# these are harmonizable but not harmonized yet!! 
harmonizable.vars <- c(
  "SURVEY_AGE", "FEMALE", 
  "USABORN_REV", "USABORNMOTHER_REV", "USABORNFATHER_REV", 
  "ETHNICITY_REV", "EDUCATION_REV", 
  "EMPLOYMENT_FULL_TIME_EMPLOYED", "EMPLOYMENT_PART_TIME_EMPLOYED", 
  "EMPLOYMENT_RETIRED", 
  "MARITALSTATUS", "INCOME", "SIZEOFHH", "INCOME_PP", 
  "GENERALHEALTH", "SR_BMI", "SR_DIABETES", "SR_HYP", "SR_CHF", 
  "SMOKING_STATUS")
# removed ALCOHOL_BINGE

harm_missingsummary <- data.frame(varname = harmonizable.vars)
row.names(harm_missingsummary) <- harmonizable.vars
for (i in harmonizable.vars) {
  harm_missingsummary[i, "pctmiss"] <- 
    100 * sum(is.na(pre_MI_data[, i])) / nrow(pre_MI_data)
  # print(i)
  # print(table(RPGEH[, i], exclude = NULL))
}

harm_missingsummary <- 
  harm_missingsummary[order(harm_missingsummary$pctmiss), ]
harm_missingsummary

## auxiliary variables ----
aux.vars <- c(
  "SURVEY_LANGUAGE", 
  "SR_TOTAL_HEIGHT_M", "SR_WEIGHT_KG", 
  "SR_STROKE", "SR_ANGINA", "SR_CANCER1", "SR_DEM", "SR_DEPRESS"
)

aux_missingsummary <- data.frame(varname = aux.vars)
row.names(aux_missingsummary) <- aux.vars
for (i in aux.vars) {
  aux_missingsummary[i, "pctmiss"] <- 
    100 * sum(is.na(pre_MI_data[, i])) / nrow(pre_MI_data)
  # print(i)
  # print(table(RPGEH[, i], exclude = NULL))
}

aux_missingsummary <- aux_missingsummary[order(aux_missingsummary$pctmiss), ]
aux_missingsummary

impute.vars <- c(harmonizable.vars, aux.vars)
# sum(duplicated(impute.vars))

# set class for variables to impute ----
# pre_MI_data_ID <- pre_MI_data %>% pull(SUBJID)
pre_MI_data <- pre_MI_data %>% select(SUBJID, all_of(impute.vars))

glimpse(pre_MI_data)

cont.vars <- c(
  "SURVEY_AGE", "INCOME_PP", 
  "SR_BMI", "SR_TOTAL_HEIGHT_M", "SR_WEIGHT_KG"
)
binary.vars <- c(
  "FEMALE", "USABORN_REV", "USABORNMOTHER_REV", "USABORNFATHER_REV", 
  "EMPLOYMENT_FULL_TIME_EMPLOYED", "EMPLOYMENT_PART_TIME_EMPLOYED", 
  "EMPLOYMENT_RETIRED", 
  "SR_DIABETES", "SR_HYP", "SR_CHF", 
  "SR_STROKE", "SR_ANGINA", "SR_CANCER1", "SR_DEM", "SR_DEPRESS"
)
ordinal.vars <- c("EDUCATION_REV", "INCOME", "SIZEOFHH", "GENERALHEALTH")
# removed alcohol binge
cat.vars <- c("ETHNICITY_REV", "MARITALSTATUS", "SMOKING_STATUS", "SURVEY_LANGUAGE")

# check
length(unique(c(cont.vars, binary.vars, ordinal.vars, cat.vars))) == length(impute.vars)

pre_MI_data <- pre_MI_data %>% 
  mutate(across(all_of(binary.vars), factor)) %>%
  mutate(across(all_of(ordinal.vars), ~factor(.x, ordered = T))) %>%
  mutate(across(all_of(cat.vars), ~factor(.x, ordered = F)))

# saveRDS(pre_MI_data, 
#         file = paste0(path_to_box, path_to_datasets, 
#                       "imputed_data/pre_MI_data_08312022.RDS"))
# check
glimpse(pre_MI_data)

# initiate imputation ----

# use proper methods
ini <- mice(pre_MI_data, maxit = 0, 
            defaultMethod = c("pmm", "logreg", "polyreg", "polr"), seed = 12345)
meth <- ini$method
pred <- ini$predictorMatrix

# keep SUBJID in the dataset but don't use it for imputation
pred[, "SUBJID"] <- 0
pred["SUBJID", ] <- 0
# change predictor matrix so income_pp doesn't predict hhincome or sizeofhh
pred[c("INCOME", "SIZEOFHH"), "INCOME_PP"] <- 0

# # change predictor matrix so H_usborn does not predict the following variables
# # this is based on the the logged events of the initial imputation
# # MICE threw out H_usborn when imputing these variables 
# pred[c("H_edu_4", "H_marit", "H_hhincome_cat", "H_hhsize_cat", "H_income_pp",       
#        "H_health_5", "H_bmi", "H_alcbinge_monthly", "H_smk", "SR_TOTAL_HEIGHT_M", 
#        "SR_WEIGHT_KG"), "H_usborn"] <- 0

# run imputation and diagnostics ---- 
imp_fcs <- mice(pre_MI_data, m = 40, maxit = 10, pred = pred, meth = meth, 
                defaultMethod = c("pmm", "logreg", "polyreg", "polr"), 
                seed = 12345)
save(imp_fcs,
     file = paste0(path_to_box, path_to_datasets,
                   "imputed_data/imp_fcs_08312022.R"))

# load the mice object post MI
# load(file = paste0(path_to_box, path_to_datasets,
#                    "imputed_data/imp_fcs_08312022.R"))

# examine diagnostics
imp_fcs$loggedEvents
plot(imp_fcs)

# prep and save stacked imputed dataset ----
imp.temp <- list()
for (i in 1:imp_fcs$m){
  imp.temp[[i]] <- complete(imp_fcs, action = i) %>% 
    mutate(imp = i)
}

imp_data_stacked <- do.call(rbind, imp.temp)

# saveRDS(imp_data_stacked,
#         file = paste0(path_to_box, path_to_datasets,
#                       "imputed_data/imp_fcs_08312022_stacked.RDS"))
