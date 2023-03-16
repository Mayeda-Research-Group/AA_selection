# variable harmonization between CHIS and RPGEH
# restrict to (1) subjects 60+ and top code at 85, (2) insured only, 
# and (3) Asians only

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "table1", "labelled"
       # "mice", "rlang"
       # "magrittr", "foreign", "ggplot2",
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "openxlsx", "ggpubr", "mitools", "lmtest"
)

options(scipen = 999, digits = 8)

# set up paths and data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_proj <- "Asian_Americans_dementia/Manuscripts/AA_selection/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
rpgeh_imp <- readRDS(paste0(path_to_box, path_to_datasets, 
                            "imputed_data/imp_fcs_08312022_stacked.RDS"))

# variable harmonization ----

glimpse(rpgeh_imp)

rpgeh_imp <- rpgeh_imp %>% 
  mutate(
    H_age = ifelse(SURVEY_AGE > 85, 85, SURVEY_AGE), # top code at 85
    H_female = ifelse(FEMALE == 1, 1, 0), 
    H_usborn = ifelse(USABORN_REV == 1, 1, 0),
    H_usborn_m = ifelse(USABORNMOTHER_REV == 1, 1, 0),
    H_usborn_f = ifelse(USABORNFATHER_REV == 1, 1, 0),
    # H_gen = case_when(
    #   H_usborn == 0 ~ 1, 
    #   H_usborn == 1 & (H_usborn_m == 0 | H_usborn_f == 0) ~ 2, 
    #   H_usborn == 1 & H_usborn_m == 1 & H_usborn_f == 1 ~ 3
    # ),
    H_ethn = case_when(
      ETHNICITY_REV %in% c(1:8) ~ as.numeric(ETHNICITY_REV),
      ETHNICITY_REV == 9 ~ 10, 
      ETHNICITY_REV == 10 ~ 9,
      TRUE ~ NA_real_
    ),
    H_edu_4 = case_when(
      EDUCATION_REV %in% c(1,2) ~ 1,
      EDUCATION_REV == 3 ~ 2,
      EDUCATION_REV == 4 ~ 3, 
      EDUCATION_REV %in% c(5,6) ~ 4
    ),  
    H_work = ifelse(
      EMPLOYMENT_FULL_TIME_EMPLOYED == 1 | EMPLOYMENT_PART_TIME_EMPLOYED == 1, 
      1, 0), 
    H_retired = ifelse(EMPLOYMENT_RETIRED == 1, 1, 0),
    H_marit = ifelse(MARITALSTATUS == 2, 1, 0),
    hhincome_top = case_when(
      INCOME == 1 ~ 20000,
      INCOME == 2 ~ 39999, 
      INCOME == 3 ~ 59999, 
      INCOME == 4 ~ 99999, 
      INCOME == 5 ~ 139999,
      TRUE ~ NA_real_
    ), 
    H_hhincome_5 = as.numeric(INCOME), 
    H_hhsize_5 = as.numeric(SIZEOFHH),
    H_hhsize_3 = ifelse(H_hhsize_5 >= 3, 3, H_hhsize_5), 
    H_income_pp = hhincome_top / sqrt(H_hhsize_5), 
    H_health_5 = as.numeric(GENERALHEALTH),
    H_health_3 = case_when(
      GENERALHEALTH %in% c(1,2) ~ 1,
      GENERALHEALTH == 3 ~ 2,
      GENERALHEALTH %in% c(4,5) ~ 3
    ), 
    H_bmi = SR_BMI,
    H_diab = ifelse(SR_DIABETES == 1, 1, 0), 
    H_hyp = ifelse(SR_HYP == 1, 1, 0),
    H_smk = as.numeric(SMOKING_STATUS), 
    # add sample weight
    smplwt = 1, 
    # add study indicator
    in_rpgeh = 1, 
    in_chis = 0
  )

glimpse(rpgeh_imp)

# keep only the harmonized variables 
rpgeh_imp <- rpgeh_imp %>% 
  select(SUBJID, starts_with("H_"), smplwt, in_rpgeh, in_chis, imp) %>% 
  rename(ID = SUBJID) %>% 
  mutate(ID = as.character(ID))

# save ----
# saveRDS(rpgeh_imp, 
#         paste0(path_to_box, path_to_datasets, "rpgeh_imp_harmonized.RDS"))
