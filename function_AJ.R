estimate_AJ <- function(data) {
  # data <- boot_data
  require(survival)
  require(broom)
  
  # because of bootstrap resampling, the original ID's are no longer unique
  # need to create new ID's for survfit
  data <- data %>% mutate(new_ID = 1:n())

  crude_AJ <- survfit(Surv(SURVEY_AGE, MAIN_DEM_V1_END_AGE, endtype) ~ 1, 
                      id = new_ID, data = data) %>% 
    tidy() %>% 
    filter(state == "dementia") %>% 
    mutate(age_grp = cut(time, 
                         breaks = c(60, 65, 70, 75, 80, 85, 90, 120), 
                         labels = c(60, 65, 70, 75, 80, 85, 90),
                         right = FALSE)) %>% 
    group_by(age_grp) %>% 
    slice_min(n = 1, order_by = time) %>% 
    ungroup() %>% 
    filter(age_grp != 60) %>%
    select(time, estimate) 
  
  wt_AJ <- survfit(Surv(SURVEY_AGE, MAIN_DEM_V1_END_AGE, endtype) ~ 1, 
                   id = new_ID, weights = sw1, data = data) %>% 
    tidy() %>% 
    filter(state == "dementia") %>% 
    mutate(age_grp = cut(time, 
                         breaks = c(60, 65, 70, 75, 80, 85, 90, 120), 
                         labels = c(60, 65, 70, 75, 80, 85, 90),
                         right = FALSE)) %>% 
    group_by(age_grp) %>% 
    slice_min(n = 1, order_by = time) %>% 
    ungroup() %>% 
    filter(age_grp != 60) %>%
    select(time, estimate) 
  
  out <- c(crude_AJ$estimate, round(crude_AJ$time, 2), 
           wt_AJ$estimate, round(wt_AJ$time, 2))
  
  names(out) <- c(
    paste0("crude_AJ_est_", seq(65, 90, 5)), 
    paste0("crude_AJ_time_", seq(65, 90, 5)),
    paste0("wt_AJ_est_", seq(65, 90, 5)), 
    paste0("wt_AJ_time_", seq(65, 90, 5))
  )
  
  return(out)
}
  
# function testing ----

# library(tidyverse)
# library(boot)
# library(survival)
# library(broom)
# 
# path_to_box <- "~/Library/CloudStorage/Box-Box/"
# path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
# load(paste0(path_to_box, path_to_datasets, "bootstrap_prep_data.RData"))
# 
# # carry out analysis as an example 
# boot_data <- chis_rpgeh %>% filter(H_ethn == "Japanese", imp %in% c(0, 1))
# p_rpgeh <- sum(boot_data$in_rpgeh) / sum(boot_data$smplwt)
# uncon_rpgeh_odds <- p_rpgeh / (1 - p_rpgeh)
# mod <- glm(fmls["Japanese"], data = boot_data,
#            family = binomial(link = "logit"), weights = boot_data$smplwt)
# boot_data <- boot_data %>% 
#   mutate(
#     # calculate predicted probabilities and weights
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
# # now call the estimate_AJ function on boot_data
# estimate_AJ(boot_data)
# # Warining messages are usually due to estimating CI limits in survfit,
# # but we only extract the point estimates, so it's not a problem. 

