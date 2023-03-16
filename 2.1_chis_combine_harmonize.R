# combine 2005, 2007, and 2009 CHIS data

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
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_proj <- "Asian_Americans_dementia/Manuscripts/AA_selection/"

chis_2005 <- read_sas(
  paste0(path_to_box, path_to_proj, 
         "CHIS_2005_2009_data/CHIS_2005_ADULT/adult.sas7bdat")
)
chis_2007 <- read_sas(
  paste0(path_to_box, path_to_proj, "CHIS_2007_data/Adult SAS/adult.sas7bdat")
)
chis_2009 <- read_sas(
  paste0(path_to_box, path_to_proj, 
         "CHIS_2005_2009_data/CHIS_2009_ADULT/adult.sas7bdat")
)

# add restrictions ----
# insured, age 60+, and Asian
# select variables used in harmonization and weight generation
vars <- c(
  "PUF_ID", "SRAGE_P", "SRSEX", 
  "AH33NEW", "AH34NEW", "AH35NEW", # nativity
  "AHEDUC", "WRKST", "AK2", # education, work, retired
  "MARIT2", 
  "AK22_P", "HHSIZE_P", # household income and size
  "AB1", "AB22", "AB29", # general health, diabetes, hypertension
  "BMI_P", "SMOKING", 
  "RAKEDW0",
  "INTVLANG"
)

chis_res_2005 <- chis_2005 %>% filter(SRAGE_P >= 60, INS == 1) %>% 
  mutate(
    H_ethn = case_when(
      SRPI == 1 & ASIAN9 == -1 ~ 8, # PI
      SRPI == 1 & ASIAN9 != -1 ~ 9, # multiple
      ASIAN9 == 5 ~ 1, # South Asian
      ASIAN9 == 1 ~ 2, # Chinese
      ASIAN9 == 2 ~ 3, # Japanese
      ASIAN9 == 3 ~ 4, # Korean
      ASIAN9 == 4 ~ 5, # Filipino
      ASIAN9 == 6 ~ 6, # Vietnamese
      ASIAN9 %in% c(7,8) ~ 7, # Cambodian + Other SE Asian
      ASIAN9 == 9 ~ 9, # multiple
      TRUE ~ 10 # all else
    )
  ) %>% 
  filter(H_ethn != 10) %>% # restrict to Asians only
  select(H_ethn, all_of(vars))

chis_res_2007 <- chis_2007 %>% filter(SRAGE_P >= 60, INS == 1) %>% 
  mutate(
    H_ethn = case_when(
      SRPI == 1 & ASIAN9 == -1 ~ 8, # PI
      SRPI == 1 & ASIAN9 != -1 ~ 9, # multiple
      ASIAN9 == 5 ~ 1, # South Asian
      ASIAN9 == 1 ~ 2, # Chinese
      ASIAN9 == 2 ~ 3, # Japanese
      ASIAN9 == 3 ~ 4, # Korean
      ASIAN9 == 4 ~ 5, # Filipino
      ASIAN9 == 6 ~ 6, # Vietnamese
      ASIAN9 %in% c(7,8) ~ 7, # Other SE Asian
      ASIAN9 == 9 ~ 9, # multiple
      TRUE ~ 10 # all else
    )
  ) %>% 
  filter(H_ethn != 10) %>% 
  select(H_ethn, all_of(vars))

chis_res_2009 <- chis_2009 %>% filter(SRAGE_P >= 60, INS == 1) %>% 
  mutate(
    H_ethn = case_when(
      SRPI == 1 & ASIAN8 == -1 ~ 8, # PI
      SRPI == 1 & ASIAN8 != -1 ~ 9, # multiple
      ASIAN8 == 5 ~ 1, # South Asian
      ASIAN8 == 1 ~ 2, # Chinese
      ASIAN8 == 2 ~ 3, # Japanese
      ASIAN8 == 3 ~ 4, # Korean
      ASIAN8 == 4 ~ 5, # Filipino
      ASIAN8 == 6 ~ 6, # Vietnamese
      ASIAN8 == 7 ~ 7, # Cambodian + Other SE Asian
      ASIAN8 == 8 ~ 9, # multiple
      TRUE ~ 10 # all else
    )
  ) %>% 
  filter(H_ethn != 10) %>% # restrict to Asians only
  select(H_ethn, all_of(vars))

# combine datasets ----
chis_combine <- rbind(chis_res_2005, chis_res_2007, chis_res_2009)

# harmonize variables ----

chis_combine <- chis_combine %>% 
  mutate(
    H_age = SRAGE_P, 
    H_female = ifelse(SRSEX == 2, 1, 0), 
    H_usborn = ifelse(AH33NEW == 1, 1, 0), 
    H_usborn_m = case_when(
      AH34NEW == 1 ~ 1,
      AH34NEW == 2 ~ 0,
      TRUE ~ NA_real_
    ), 
    H_usborn_f = case_when(
      AH35NEW == 1 ~ 1,
      AH35NEW == 2 ~ 0,
      TRUE ~ NA_real_
    ), 
    # H_gen = case_when(
    #   H_usborn == 0 ~ 1, 
    #   H_usborn == 1 & (H_usborn_m == 0 | H_usborn_f == 0) ~ 2, 
    #   H_usborn == 1 & H_usborn_m == 1 & H_usborn_f == 1 ~ 3, 
    #   TRUE ~ NA_real_
    # ),
    H_edu_4 = case_when(
      AHEDUC %in% c(91, 1, 2) ~ 1, 
      AHEDUC == 3 ~ 2,
      AHEDUC %in% 4:5 ~ 3,
      AHEDUC %in% 6:10 ~ 4
    ),
    H_work = ifelse(WRKST %in% c(1,2,3), 1, 0),
    H_retired = ifelse(AK2 == 5, 1, 0),
    H_marit = ifelse(MARIT2 %in% c(1,2), 1, 0),
    hhincome_top = case_when(
      AK22_P < 20000 ~ 20000,
      AK22_P < 40000 ~ 39999, 
      AK22_P < 60000 ~ 59999, 
      AK22_P < 100000 ~ 99999, 
      AK22_P >= 100000 ~ 139999, 
      TRUE ~ NA_real_
    ), 
    H_hhincome_5 = case_when(
      AK22_P < 20000 ~ 1,
      AK22_P < 40000 ~ 2, 
      AK22_P < 60000 ~ 3, 
      AK22_P < 100000 ~ 4, 
      AK22_P >= 100000 ~ 5, 
      TRUE ~ NA_real_
    ), 
    H_hhsize_5 = ifelse(HHSIZE_P >= 5, 5, HHSIZE_P), # used for income_pp
    H_hhsize_3 = ifelse(H_hhsize_5 >= 3, 3, H_hhsize_5), # used in model generation
    H_income_pp = hhincome_top / sqrt(H_hhsize_5),
    H_health_3 = case_when(
      AB1 %in% c(1,2) ~ 1,
      AB1 == 3 ~ 2,
      AB1 %in% c(4,5) ~ 3
    ),
    H_health_5 = AB1,
    H_bmi = BMI_P,
    H_diab = ifelse(AB22 %in% c(1,3), 1, 0),
    H_hyp = ifelse(AB29 == 1, 1, 0), 
    H_smk = case_when(
      SMOKING == 3 ~ 1,
      SMOKING == 2 ~ 2, 
      SMOKING == 1 ~ 3
    ), 
    # adjust weights for pooling 3 years of CHIS
    smplwt = RAKEDW0/3,
    # add study indicators
    in_chis = 1, 
    in_rpgeh = 0
  ) 

glimpse(chis_combine)

# keep only the harmonized variables 
chis_combine <- chis_combine %>% 
  select(PUF_ID, starts_with("H_"), smplwt, in_chis, in_rpgeh, INTVLANG) %>% 
  rename(ID = PUF_ID) 

# save ----
# saveRDS(chis_combine, 
#         paste0(path_to_box,"Asian_Americans_dementia_data/aa_selection/", 
#                "chis_2005_to_2009_harmonized.RDS"))

