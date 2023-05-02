# bootstrap script to calculate incidence rates by ethnicity
# including unweighted and weighted crude and age-adjusted IR's

# Hoffman setup ----
# Read in the arguments listed in the:
# R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  
## expression:
args = (commandArgs(TRUE))

# Check to see if arguments are passed and set default values if not.
# If so, parse the arguments. (I only have one argument here.)
if (length(args) == 0) {
  print("No arguments supplied.")
  ##supply default values
  scenario_num <- 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
# print values just to make sure:
print(scenario_num)

# also make sure that scenario_num is not greater than 4
if (scenario_num > 4) {
  print("out of bound senario number")
  # supply default value
  scenario_num <- 1
}

# we will use scenario_num to determine 
# (1) which imputed datasets we loop the bootstrapping process over
# e.g. senario_num = 1, then imp = 1:10
# senario_num = 2, then imp = 11:20, etc. so that we have 4 jobs in the array
# (2) which seed we start out with before the loop

# set up packages ----
## for desktop ----
# if (!require("pacman")) 
#   install.packages("pacman", repos='http://cran.us.r-project.org')
# 
# p_load("haven", "tidyverse", "twang", "boot"
#        # "mice", rlang",
#        # "magrittr", "foreign", "ggplot2", "table1", "labelled"
#        # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
#        # "openxlsx", "lmtest", "mitools", "ggpubr", "patchwork", "RColorBrewer"
# )
# 
# options(scipen = 999, digits = 8)

## for Hoffman ----
library(tidyverse)
library(boot)

# set up paths and data ----

### when test running on desktop
# path_to_box <- "~/Library/CloudStorage/Box-Box/"
# path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
# path_to_output <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/output/"

# load the prepped data generated in 4.1 and the IR calculation function
# load(paste0(path_to_box, path_to_datasets, "bootstrap_prep_data.RData"))
# source(paste0(path_to_box, 
#               "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
#               "cleaned_scripts/function_IR_calc.R"))

### when running on Hoffman
# load the prepped data generated in 4.1 and the IR calculation function
load("/u/home/y/yixzhou/AA_selection/bootstrap_prep_data.RData")
source("/u/home/y/yixzhou/AA_selection/function_IR_calc.R")

# overall description: 
# take one imputed copy of RPGEH with CHIS
# iterate over the ethnicities by 
# isolating one ethnicity and running bootstrap sampling 1100 times
# collect (1) IR output into csv output files by ethnicity
# and (2) warnings and running time by ethnicity and imputation

# bootstrapping function ----

# ir_boot calculates one line of statistics for one ethnicity and one imputation
# (1) unweighted and weighted age specific IR's - 12 statistics
# (2) unweighted and weighted pys and cases - 4 statistics
# (3) unweighted and weighted crude and age-adjusted IR's - 4 statistics
# (4) auxiliary variables for bootstrapping: 
#     t = index of bootstrap sample, imp = index of imputed dataset used, 
#     warning = indicator for potential glm warning 

ir_boot <- function(data, formula, indices) {
  # input data should be filtered to a certain ethnicity
  boot_data <- data[indices,]

  # calculate unconditional RPGEH sampling odds
  p_rpgeh <- sum(boot_data$in_rpgeh) / sum(boot_data$smplwt)
  uncon_rpgeh_odds <- p_rpgeh / (1 - p_rpgeh)

  # fit the model
  # tryCatch records if there is a glm warning and 
  # accumulates the number of warnings
  tryCatch(
    {
      glm(as.formula(formula), data = boot_data,
          family = binomial(link = "logit"), weights = boot_data$smplwt)
      ww <<- c(ww, 0)
    },
    warning = function(w){
      n_w <<- n_w + 1
      ww <<- c(ww, 1)
    })
  
  mod <- glm(as.formula(formula), data = boot_data,
             family = binomial(link = "logit"), weights = boot_data$smplwt)

  boot_data <- boot_data %>% 
    mutate(
      # calculate predicted probabilities and weights
      p1 = predict.glm(mod, type = "response"), 
      iow1 = (1 - p1) / p1,
      sw1 = ifelse(in_rpgeh == 1, iow1 * uncon_rpgeh_odds, 1),
      # sw1_twang = sw1 * smplwt # this is only for twang - bal.stat
    ) %>% 
    # we don't need CHIS data for IR calculations
    filter(in_rpgeh == 1) %>% 
    # join RPGEH pys data 
    left_join(., pys_data, "ID")
  
  # calculate IR to output
  out <- calculate_ir(
    data = boot_data, 
    age_cat = age_cat_labels, 
    # this is the population we standardize to 
    std_pop = c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
  )
  
  return(out)
}
  

# bootstrap starts here ----
ethns_list <- names(fmls)
# i <- 1 # which imputation being worked on, will be an argument

# set up the index of imputations using scenario_num
# scenario_num <- 4
imp_start <- (scenario_num - 1) * 10 + 1
imp_end <- scenario_num * 10
seed <- scenario_num * 1234
set.seed(seed)

for (i in imp_start:imp_end) {
  for (ethn in ethns_list) {
    start <- Sys.time() # record starting time
    n_w <- 0 # a counter for glm warnings
    ww <- numeric() # an index tracker for which sample has warning
    
    # subset to the specific ethnicity and the i-th imputed RPGEH dataset
    input_data <- chis_rpgeh %>%
      filter(H_ethn == ethn, imp %in% c(0, i))
    
    # run bootstrapping
    res <- boot(
      data = input_data,
      statistic = ir_boot,
      R = 1100,
      strata = input_data$in_rpgeh, # added on Mar 23 for stratified resampling
      formula = fmls[ethn]
    )
    
    end <- Sys.time() # record ending time
    
    # save bootstrapping results to output
    # res$t0 is results using full dataset, no resampling
    # res$t is results with resampling
    out_df <- rbind(res$t0, res$t)
    # append with sampling number and imputation number (i)
    out_df <- cbind(out_df, t = c(0:res$R), imp = i, warning = ww) %>% as_tibble()
    
    write_csv(out_df, 
              append = i != 1, # only create new file when i = 1
              file = paste0(
                # desktop path
                # path_to_box, path_to_output, 
                # "bootstrap_results/bootstrap_ir_", ethn, ".csv"
                "/u/home/y/yixzhou/AA_selection/output/", 
                "bootstrap_ir_", ethn, "_w_warn_index.csv"
              ))
    
    # save the number of warnings and duration of the run for review
    w_df <- tibble(ethnicity = ethn, n_warns = n_w, imp = i, 
                   duration = difftime(end, start, units = 'mins'))
    write_csv(w_df, 
              # only create new file when it's the first ethnicity and i = 1
              append = !(ethn == ethns_list[1] & i == 1),
              file = paste0(
                # desktop path
                # path_to_box, path_to_output, 
                # "bootstrap_results/boostrap_warnings.csv"
                "/u/home/y/yixzhou/AA_selection/output/", 
               "bootstrap_warnings_w_warn_index.csv"
                ))
  }
}
