# table 1 preparation
# for CHIS and RPGEH (pre and post imputation)

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "table1", "labelled", "openxlsx", "survey", 
       "gtsummary"
       # "mice", rlang", "twang", 
       # "magrittr", "foreign", "ggplot2", 
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "lmtest", "ggpubr",  "mitools", 
)

options(scipen = 999, digits = 8)

# set up paths and data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
path_to_output <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/output/"

chis_harm <- readRDS(paste0(path_to_box, path_to_datasets, 
                            "chis_2005_to_2009_harmonized.RDS"))
rpgeh_harm_imp <- readRDS(paste0(path_to_box, path_to_datasets, 
                                 "rpgeh_imp_harmonized.RDS")) 
rpgeh_tte <- read_sas(paste0(path_to_box,
                             "Asian_Americans_dementia_data/analysis_data_tables/",
                             "aa_adrd_cardiometabolic_tte.sas7bdat"))
rpgeh_pre_MI <- readRDS(paste0(path_to_box, path_to_datasets, 
                               "imputed_data/pre_MI_data_08312022.RDS"))

source(paste0(path_to_box, 
              "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
              "cleaned_scripts/functions_summary.R"))

# table 1 for CHIS ----
# apply format
chis_harm <- chis_harm %>% t1_relabel()

# set up svydesign objects
svy.d <- chis_harm %>% 
  select(
    H_age, H_female, H_usborn, H_edu_4,
    H_marit, H_hhsize_3, H_income_pp, H_work, H_retired, 
    H_health_3, H_bmi, H_diab, H_hyp, H_smk, H_ethn, smplwt) %>% 
  svydesign(~ 1, data = ., weights = ~ smplwt)

# combined table 1 ----
# where for continuous variables, the weighted mean and SD are used
# and for categorical variables, the weighted percentages are used

chis_combined <- tbl_svysummary(
  data = svy.d, 
  by = H_ethn, 
  statistic = list(all_continuous() ~ "{mean} ({sd})", 
                   all_categorical() ~ "{p}"), 
  missing_text = "Missing", 
  digits = list(H_age ~ 1, H_bmi ~ 1, H_income_pp ~ 0, 
                all_categorical() ~ 1)
)

chis_counts_ethn <- chis_harm %>% 
  count(H_ethn) %>% 
  rename(`Unweighted N` = n) %>% 
  left_join(
    chis_harm %>% count(H_ethn, wt = smplwt) %>% 
      mutate(n = round(n)) %>% 
      rename(`CHIS survey-weighted N` = n), 
    by = "H_ethn") %>% 
  t()

write.xlsx(
  list(
    "combined" = as_tibble(chis_combined),
    "counts" = chis_counts_ethn
  ),
  rowNames = TRUE,
  file = paste0(path_to_box, path_to_output, "table1s_without_n/chis_tbl1.xlsx")
)

## overall summary statistics for age ----
svymean(~H_age, svy.d) # mean = 71.0
svyvar(~H_age, svy.d) # var = 53.6
sqrt(svyvar(~H_age, svy.d)[1]) # sd = 7.3

# table 1 for RPGEH ----

## pre-imputation ----
names(rpgeh_pre_MI)

# merge in FU time and FU end type
rpgeh_pre_MI <- rpgeh_tte %>% 
  select(SUBJID, MAIN_DEM_V1_FU_TIME, MAIN_DEM_V1_END_TYPE) %>% 
  right_join(rpgeh_pre_MI, by = "SUBJID") %>% 
  t1_relabel_pre_MI()

### overall summary statistics ----
# survey age in pre-MI dataset is not topcoded
mean(rpgeh_pre_MI$SURVEY_AGE) # 70.3
sd(rpgeh_pre_MI$SURVEY_AGE) # 7.1

mean(rpgeh_pre_MI$MAIN_DEM_V1_FU_TIME) # 9.9
sd(rpgeh_pre_MI$MAIN_DEM_V1_FU_TIME) # 4.6

tbl1_pre_MI <- rpgeh_pre_MI %>%
  table1(
    data = .,
    ~ SURVEY_AGE + FEMALE + USABORN_REV + 
      SURVEY_LANGUAGE + EDUCATION_REV + 
      EMPLOYMENT_FULL_TIME_EMPLOYED + EMPLOYMENT_PART_TIME_EMPLOYED +
      EMPLOYMENT_RETIRED +
      MARITALSTATUS + SIZEOFHH + # INCOME_PP + 
      GENERALHEALTH + SR_BMI + SR_DIABETES + SR_HYP + SMOKING_STATUS +
      MAIN_DEM_V1_FU_TIME + MAIN_DEM_V1_END_TYPE + 
      # auxiliary variables used in imputation 
      USABORNMOTHER_REV + USABORNFATHER_REV + 
      SR_TOTAL_HEIGHT_M + SR_WEIGHT_KG +
      SR_CHF + SR_STROKE + SR_ANGINA + SR_CANCER1 +
      SR_DEM + SR_DEPRESS
    | ETHNICITY_REV, 
    overall = NULL, 
    render.continous = my.render.cont, 
    render.categorical = "PCT"
    )

 
tbl1_pre_MI_incomepp <- rpgeh_pre_MI %>% #  t1_relabel_pre_MI() %>% 
  group_by(ETHNICITY_REV) %>% 
  summarize(mean = mean(INCOME_PP, na.rm = TRUE) %>% round(digits = 0), 
            SD = sd(INCOME_PP, na.rm = TRUE) %>% round(digits = 0), 
            n_missing = sum(is.na(INCOME_PP)), 
            perc_missing = round(sum(is.na(INCOME_PP) / n())*100, digits = 1))

# write.xlsx(list(tbl1_pre_MI, tbl1_pre_MI_incomepp), rowNames = TRUE,
#            file = paste0(path_to_box, path_to_output,
#                          "table1s_without_n/RPGEH_pre_MI_tbl1.xlsx"))

### additional summary on race/ethnicity ----
race_ethn_data <- rpgeh_tte %>% 
  select(SUBJID, starts_with("ETHN_")) %>% 
  right_join(rpgeh_pre_MI, by = "SUBJID") %>% 
  mutate(
    Multiracial = case_when(
      ETHN_WHITE == 0 &   
        ETHN_AFRICAN_AMERICAN == 0 & 
        ETHN_AFRICAN == 0 & 
        ETHN_AFRO_CARIBBEAN == 0 & 
        ETHN_MEXICAN == 0 & 
        ETHN_CENTRAL_S_AMERICAN == 0 & 
        ETHN_PUERO_RICAN == 0 & 
        ETHN_CUBAN == 0 & 
        ETHN_LATINO == 0 & 
        ETHN_MIDDLE_EASTERN == 0 & 
        ETHN_ASHKENAZI_JEWISH == 0 & 
        ETHN_OTHER == 0 ~ 0, # "No non-Asian ethnicity", 
      TRUE ~ 1 # "At last one non-Asian ethnicity"
    ) %>% 
      as.logical(),
    across(starts_with("ETHN_"), as.logical)
  )

table1(~ Multiracial | ETHNICITY_REV, 
       data = race_ethn_data, 
       render = logical.rndr, render.categorical = "PCT") 
  # write.xlsx(rowNames = TRUE,
  #            file = paste0(path_to_box, path_to_output,
  #                          "table1s_without_n/RPGEH_multiracial.xlsx"))

race_ethn_data %>% 
  filter(ETHNICITY_REV == "Multiple Ethnicities") %>%
  table1(~ ETHN_SOUTH_ASIAN + ETHN_CHINESE + ETHN_JAPANESE + ETHN_KOREAN + 
         ETHN_FILIPINO + ETHN_VIETNAMESE + ETHN_ANY_PAC_ISLANDER + 
         ETHN_OTHER_SE_ASIAN | ETHNICITY_REV, 
       data = ., overall = FALSE, render = logical.rndr)


## post-imputation ----
rpgeh_harm_imp <- rpgeh_harm_imp %>% 
  t1_relabel() %>% 
  mutate(SUBJID = as.double(ID)) %>%
  # pull SURVEY_LANGUAGE from full dataset
  left_join(., select(rpgeh_tte, SUBJID, SURVEY_LANGUAGE), by = "SUBJID") %>% 
  # pull FU time and end type 
  left_join(., select(rpgeh_pre_MI, SUBJID,
                      MAIN_DEM_V1_FU_TIME, MAIN_DEM_V1_END_TYPE), by = "SUBJID")
var_label(rpgeh_harm_imp) <- list(SURVEY_LANGUAGE = "Survey language")

###  variables with no missingness ----
# age, female, ethn, work, retired, diab, hyp, chf
rpgeh_nomiss <- rpgeh_harm_imp %>% filter(imp == 1) %>% 
  table1(
    ~ H_age + H_female + SURVEY_LANGUAGE + H_work + H_retired + 
      H_diab + H_hyp + 
      MAIN_DEM_V1_FU_TIME + MAIN_DEM_V1_END_TYPE | H_ethn,
    data = ., overall = NULL, render.continous = my.render.cont, 
    render.categorical = "PCT"
  )

# write.xlsx(rpgeh_nomiss, rowNames = TRUE,
#            file = paste0(path_to_box, path_to_output,
#                          "table1s_without_n/RPGEH_nomiss_tbl1.xlsx"))

### imputed categorical vars ----
# marit, usborn, usborn_m/f, gen, hhsize_3, edu_4, health_3, smk, alcbinge_monthly

rpgeh_miss_cat <- rpgeh_harm_imp %>% 
  table1(
    ~ H_marit + H_usborn + #H_gen 
      + H_hhsize_3 + H_edu_4 + H_health_3 + H_smk 
    | H_ethn,
    data = ., overall = NULL, 
    render.categorical = "PCT" # my.render.cat_imp
  )

# write.xlsx(rpgeh_miss_cat, rowNames = TRUE,
#            file = paste0(path_to_box, path_to_output,
#                          "table1s_without_n/RPGEH_imp_cat_tbl1.xlsx"))

### imputed continuous vars ----
# bmi, income_pp

rpgeh_miss_cont <- imp_cts_var_tbl1(rpgeh_harm_imp,
                                    list("H_income_pp" = 0, "H_bmi" = 1))
rownames(rpgeh_miss_cont)[2:3] <- 
  var_label(rpgeh_harm_imp)[rownames(rpgeh_miss_cont)[2:3]] %>% unlist()
rpgeh_miss_cont <- data.frame(rpgeh_miss_cont)

# write.xlsx(rpgeh_miss_cont, rowNames = TRUE,
#            file = paste0(path_to_box, path_to_output,
#                          "table1s_without_n/RPGEH_imp_cts_tbl1.xlsx"))

