# functions that set labels and formats, etc, for datasets and results

# function to re-factor datasets
# and give meaningful labels for variables

# use this for harmonized variables
t1_relabel <- function(df) {
  require(table1)
  require(labelled)
  tb1_df <- df %>% 
    mutate(
      H_female = factor(H_female, levels = c(1, 0), labels = c("Female", "Male")),
      H_ethn = factor(
        H_ethn, levels = c(2,5,3,4,8,1,6,7,9),
        labels = c("Chinese", "Filipino", "Japanese", "Korean",
                   "Pacific Islander", "South Asian", "Vietnamese",
                   "Other SE Asian", "Multiple Ethnicities")
        # levels = 1:9, 
        # labels = c("South Asian", "Chinese", "Japanese", "Korean", 
        #            "Filipino", "Vietnamese", "Other Southeast Asian", 
        #            "Pacific Islanders", "Multiple Asian ethnicities")
                      ), 
      H_usborn = factor(H_usborn, levels = c(1, 0), labels = c("Yes", "No")),
      H_usborn_m = factor(H_usborn_m, levels = c(1, 0), labels = c("Yes", "No")),
      H_usborn_f = factor(H_usborn_f, levels = c(1, 0), labels = c("Yes", "No")),
      # H_gen = factor(
      #   H_gen, levels = 1:3, 
      #   labels = c("1st generation", "2nd generation", "3rd generation or higher")), 
      # edu_7 = factor(
      #   edu_7, levels = c(0:6),
      #   labels = c("no formal education",
      #              "Grade 1-8", "Grade 9-11", "Grade 12/HS diploma/GED",
      #              "technical/trade/vocational school/some college",
      #              "college", "graduate school")),
      H_edu_4 = factor(
        H_edu_4, levels = 1:4, 
        labels = c("<= Grade 11", "Grade 12/HS diploma/GED", 
                   "technical/trade/vocational school/some college", 
                   "college and above")), 
      H_marit = factor(H_marit, levels = c(1, 0), labels = c("Yes", "No")), 
      H_work = factor(H_work, levels = c(1, 0), labels = c("Yes", "No")),
      H_retired = factor(H_retired, levels = c(1, 0), labels = c("Yes", "No")),
      H_hhsize_3 = factor(
        H_hhsize_3, levels = 1:3, 
        labels = c("Live alone", "Two", "Three and more")),
      H_hhsize_5 = factor(
        H_hhsize_5, levels = 1:5, 
        labels = c("One", "Two", "Three", "Four", "Five+")),
      H_health_3 = factor(
        H_health_3, levels = 1:3, 
        labels = c("Excellent or very good", "Good", "Fair or poor")),
      H_health_5 = factor(
        H_health_3, levels = 1:5, 
        labels = c("Excellent", "Very good", "Good", "Fair", "Poor")),
      H_diab = factor(H_diab, levels = c(1, 0), labels = c("Yes", "No")),
      H_hyp = factor(H_hyp, levels = c(1, 0), labels = c("Yes", "No")),
      # H_chf = factor(H_chf, levels = c(1, 0), labels = c("Yes", "No")),
      # H_alcbinge_monthly = factor(H_alcbinge_monthly, levels = c(1, 0), 
      #                             labels = c("Yes", "No")),
      H_smk = factor(
        H_smk, levels = c(1:3), labels = c("Never", "Former", "Current"))
      # in_rpgeh = factor(in_rpgeh, levels = c(1, 0), labels = c("RPGEH", "CHIS")),
      # in_chis = factor(in_rpgeh, levels = c(1, 0), labels = c("RPGEH", "CHIS")),
    )
    
  var_label(tb1_df) <- list(
    H_age = "Age, years [mean (SD)]", 
    H_female = "Sex/gender, n (%)", 
    H_ethn = "Asian ethnicity", 
    H_usborn = "US born, n (%)", 
    H_usborn_m = "Mother US born, n (%)", 
    H_usborn_f = "Father US born, n (%)", 
    # H_gen = "Generation, n (%)", 
    # H_edu_7 = "Education attainment", 
    H_edu_4 = "Education attainment, n (%)", 
    H_marit = "Married or living as married, n (%)", 
    H_work = "Currently working, n (%)", 
    H_retired = "Currently retired, n (%)", 
    H_hhsize_5 = "Household size (5 levels), n (%)", 
    H_hhsize_3 = "Household size (3 levels), n (%)", 
    H_income_pp = "Household adjusted income, dollars [mean (SD)]", 
    H_health_3 = "General health (3 levels), n (%)", 
    H_health_5 = "General health (5 levels), n (%)", 
    H_bmi = "BMI, mean (SD)", 
    H_diab = "Ever had diabetes, n (%)", 
    H_hyp = "Ever had hypertension, n (%)", 
    # H_chf = "Ever had congestive heart failure", 
    # H_alcbinge_monthly = "Alcohol binge monthly or more", 
    H_smk = "Smoking status, n (%)"
  )
  return(tb1_df)
}

# use this for pre-MI original RPGEH dataset

t1_relabel_pre_MI <- function(df) {
  require(table1)
  require(labelled)
  tb1_df <- df %>% 
    mutate(
      FEMALE = factor(FEMALE, levels = c(1, 0), labels = c("Female", "Male")),
      ETHNICITY_REV = factor(
        ETHNICITY_REV, levels = c(2,5,3,4,8,1,6,7,10),
        labels = c("Chinese", "Filipino", "Japanese", "Korean",
                   "Pacific Islander", "South Asian", "Vietnamese",
                   "Other SE Asian", "Multiple Ethnicities")
      ), 
      USABORN_REV = factor(USABORN_REV, levels = c(1, 0), labels = c("Yes", "No")),
      USABORNMOTHER_REV = factor(
        USABORNMOTHER_REV, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      USABORNFATHER_REV = factor(
        USABORNFATHER_REV, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      EDUCATION_REV = factor(
        EDUCATION_REV, levels = c(1:6),
        labels = c("Grade school (1-8)", "Some high school (9-11)", 
                   "High school of GED",
                   "Technical/trade/vocational school/some college",
                   "College", "Graduate school")),
      MARITALSTATUS = factor(
        MARITALSTATUS, levels = 1:4, 
        labels = c("Never married", "Married or living as married", 
                   "Separated/divorced", "Widowed")), 
      EMPLOYMENT_FULL_TIME_EMPLOYED = factor(
        EMPLOYMENT_FULL_TIME_EMPLOYED, levels = c(1, 0), labels = c("Yes", "No")),
      EMPLOYMENT_PART_TIME_EMPLOYED = factor(
        EMPLOYMENT_PART_TIME_EMPLOYED, levels = c(1, 0), labels = c("Yes", "No")),
      EMPLOYMENT_RETIRED = factor(
        EMPLOYMENT_RETIRED, levels = c(1, 0), labels = c("Yes", "No")),
      SIZEOFHH = factor(
        SIZEOFHH, levels = 1:5, 
        labels = c("One", "Two", "Three", "Four", "Five or more")),
      GENERALHEALTH = factor(
        GENERALHEALTH, levels = 1:5, 
        labels = c("Excellent", "Very good", "Good", "Fair", "Poor")),
      SR_DIABETES = factor(SR_DIABETES, levels = c(1, 0), labels = c("Yes", "No")),
      SR_HYP = factor(SR_HYP, levels = c(1, 0), labels = c("Yes", "No")),
      SR_CHF = factor(SR_CHF, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SR_STROKE = factor(SR_STROKE, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SR_ANGINA = factor(SR_ANGINA, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SR_CANCER1 = factor(SR_CANCER1, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SR_DEM = factor(SR_DEM, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SR_DEPRESS = factor(SR_DEPRESS, levels = c(1, 0), labels = c("Yes", "No")), # auxiliary
      SMOKING_STATUS = factor(
        SMOKING_STATUS, levels = 1:3, labels = c("Never", "Former", "Current")), 
      MAIN_DEM_V1_END_TYPE = factor(
        MAIN_DEM_V1_END_TYPE, 
        levels = c("DEMENTIA", "DEATH", "ADMIN CENSORED", "END OF MEMBERSHIP", 
                   "CENSORED 90+"),
        labels = c("Dementia", "Death", "Administratively Censored", 
                   "End of Membership", "Censored 90+"))
    )
  
  var_label(tb1_df) <- list(
    SURVEY_AGE = "Age, years [mean (SD)]", 
    FEMALE = "Female, n (%)", 
    ETHNICITY_REV = "Race/Ethnicity", 
    USABORN_REV = "US born, n (%)", 
    USABORNMOTHER_REV = "Mother US born, n (%)", 
    USABORNFATHER_REV = "Father US born, n (%)", 
    EDUCATION_REV = "Education attainment, n (%)", 
    MARITALSTATUS = "Married or living as married, n (%)", 
    EMPLOYMENT_FULL_TIME_EMPLOYED = "Currently working full time, n (%)", 
    EMPLOYMENT_PART_TIME_EMPLOYED = "Currently working part time, n (%)", 
    EMPLOYMENT_RETIRED = "Currently retired, n (%)", 
    SIZEOFHH = "Household size, n (%)", 
    INCOME_PP = "Income per person, dollars [mean (SD)]", 
    GENERALHEALTH = "General health, n (%)", 
    SR_BMI = "BMI, mean(SD)", 
    SR_TOTAL_HEIGHT_M = "Height, meters [mean (SD)]", 
    SR_WEIGHT_KG = "Weight, kg [mean (SD)]", 
    SR_DIABETES = "Ever had diabetes, n (%)", 
    SR_HYP = "Ever had hypertension, n (%)", 
    SR_CHF = "Ever had congestive heart failure, n (%)",
    SR_STROKE = "Ever had strokes, n (%)",
    SR_ANGINA = "Ever had angina, n (%)",
    SR_CANCER1 = "Ever had cancers, n (%)",
    SR_DEM = "Ever had dementia, n (%)",
    SR_DEPRESS = "Ever had depression, n (%)",
    SMOKING_STATUS = "Smoking status, n (%)",
    MAIN_DEM_V1_END_TYPE = "End of follow up event, n (%)", 
    MAIN_DEM_V1_FU_TIME = "Follow up time, years [mean (SD)]"
  )
  return(tb1_df)
}

# customized function for rendering continuous variables in table 1's
my.render.cont <- function(x) {
  with(
    stats.apply.rounding(stats.default(x), digits = 1,
                         rounding.fn = round_pad),
    c("", "Mean (SD)" = sprintf("%s (%s)", MEAN, SD))
  )
}

# customized function for rendering categorical variables in table 1's post-MI
my.render.cat_imp <- function(x, value.prefix = "") {
  N <- length(x) / 40
  freqs <- table(x) / 40
  pcts <- round(freqs / N * 100, 1)
  freqs_formatted <- paste0(round(freqs, 0), " (", as.numeric(pcts), "%)")
  names(freqs_formatted) <- paste0(value.prefix, names(freqs))
  
  out <- c("", freqs_formatted)
  return(out)
}

# function to apply Rubin's rule to imputed continuous variables

imp_cts_var_tbl1 <- function(df, cts_vars) {
  # cts_vars is a list that looks like:
  # cts_vars <- list("H_income_pp" = 0, "H_bmi" = 1)
  # where the names of the items are variable names, and 
  # the numbers are rounding digits
  
  out <- df %>% 
    group_by(H_ethn) %>% 
    group_keys() 
  
  for (var in names(cts_vars)) {
    temp <- df %>% 
      group_by(imp, H_ethn) %>% 
      summarise(mean_imp = mean(get(var)), var_imp = var(get(var))) %>% 
      ungroup() %>% 
      group_by(H_ethn) %>% 
      summarise(mean = mean(mean_imp), 
                sd = mean(var_imp) %>% sqrt()) %>% 
      mutate(out = paste0(round_pad(mean, cts_vars[[var]]), 
                          " (", round_pad(sd, cts_vars[[var]]), ")"))
    out[var] <- temp$out
  }
  
  return(t(out))
}
