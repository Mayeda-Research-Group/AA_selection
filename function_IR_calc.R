# functions to calculate age adjusted IRs


# data being supplied should have
# ID: id variable
# fu_start: age at start of followup
# fu_end: age at end of followup regardless of event status
# event_flag: whether end of followup is due to event, 1 = event, 0 = others
# wt: optional weight variable

# age categories should be supplied as a vector
# e.g. age_cat <- c(59, 65, 70, 75, 80, 85, 120)

# standard population by age category should be supplied as a vector
# length being the same as the number of age categories of interest
# pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587)

# function that calculates person-years and case contribution by age category 
calculate_pys <- function(data, age_cat, fu_start, fu_end, event_flag) {
  # testing
  # data <- rpgeh_dem
  # age_cat <- c(59, 65, 70, 75, 80, 85, 120)
  # fu_start <- "SURVEY_AGE"
  # fu_end = "MAIN_DEM_V1_END_AGE"
  # event_flag <- "MAIN_DEM_V1_END_DEM_FLAG"
  
  age_cat_labels <- c()
  
  for (i in 1:(length(age_cat) - 1)) {
    # testing 
    # i <- 1
    
    # initialize the age interval
    cat_start <- age_cat[i]
    cat_end <- age_cat[i+1]
    # create the variable names for py and case contribution
    age_int <- paste0(cat_start, "_", cat_end)
    vars <- paste0(c("py_", "case_"), age_int)
    py_var <- vars[1]
    case_var <- vars[2]
    
    # collect the age intervals for setting up the IR output table
    age_cat_labels <- c(age_cat_labels, age_int)
    
    # calculate py and case contribution by age category
    data <- data %>% 
      mutate(
        !! py_var := case_when(
          # start of fu after age category or end of fu before age category: 
          # contribute 0 pys
          get(fu_start) > cat_end | get(fu_end) <= cat_start ~ 0, 
          
          # start of fu before age category and end of fu during age category: 
          get(fu_start) <= cat_start & get(fu_end) <= cat_end ~ 
            get(fu_end) - cat_start, 
          
          # start of fu before age category and end of fu after age category: 
          get(fu_start) <= cat_start & get(fu_end) > cat_end ~ 
            cat_end - cat_start, 
          
          # start and end of fu during age category:
          get(fu_start) <= cat_end & get(fu_end) <= cat_end ~ 
            get(fu_end) - get(fu_start), 
          
          # start of fu during age category and end of fu after age cateogory:
          get(fu_start) <= cat_end & get(fu_end) > cat_end ~ 
            cat_end - get(fu_start), 
          
          TRUE ~ NA_real_),  
        
        # if end of fu is during age category, case contribution is the dem flag. 
        # otherwise, case contribution is zero
        !! case_var := ifelse(
          get(fu_end) <= cat_end & get(fu_end) > cat_start, 
          get(event_flag), 0
        )
      )
  }
  
  return(list(data = data, age_cat_labels = age_cat_labels))
}

# function that calculates age-adjusted incidence rates 
# given data with pys and cases (generated from calculate_pys()),
# unweighted and weighted with sampling weights

calculate_ir <- function(data, age_cat_labels, std_pop) {
  # data <- boot_data
  # age_cat_labels <- age_cat_labels
  # std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587)
  
  # using 2000 US census population as the standard population
  # set up age-specific IR table
  age_spec_ir <- tibble(
    age_range = age_cat_labels, 
    pop_count = std_pop, # std pop - count by age category
    pop_total = sum(pop_count), # std pop total 
    pop_prop = pop_count / pop_total, # std pop prop by age category
    unw_pys = NA, 
    unw_cases = NA, 
    wt_pys = NA, 
    wt_cases = NA
  )
  
  for (i in 1:length(age_cat_labels)) {
    # sum up pys and cases unweighted or weighted
    # by age category
    cat <- age_cat_labels[i]
    # unweighted
    age_spec_ir[i, "unw_pys"] <- sum(data[, paste0("py_", cat)]) 
    age_spec_ir[i, "unw_cases"] <- sum(data[, paste0("case_", cat)])
    # weighted
    age_spec_ir[i, "wt_pys"] <- sum(data[, "sw1"] * data[, paste0("py_", cat)])
    age_spec_ir[i, "wt_cases"] <- sum(data[, "sw1"] * data[, paste0("case_", cat)])
  }
  
  # this multiplication step is needed to calculate the age-adjusted IR's
  # but for age-specific rates it should only be cases / pys
  # these results are later corrected in 4.4_bootstrap_processing.R
  age_spec_ir <- age_spec_ir %>% 
    mutate(
      unw_adj_ir = (unw_cases / unw_pys) * pop_prop,
      wt_adj_ir = (wt_cases / wt_pys) * pop_prop
      )
  
  # format age-specific IRs output
  age_spec_ir_out <- c(age_spec_ir$unw_adj_ir, age_spec_ir$wt_adj_ir)
  names(age_spec_ir_out) <- paste0(
    rep(c("unw_age_spec_", "wt_age_spec_"), each = length(age_cat_labels)),
    age_cat_labels)
  
  ir_out <- age_spec_ir %>% 
    summarise(across(unw_pys:wt_adj_ir, sum)) %>% 
    mutate(
      # generate crude IR, unweighted and weighted
      unw_crude_ir = unw_cases / unw_pys,
      wt_crude_ir = wt_cases/ wt_pys
      # crude_ir = cases / pys,
      # crude_ir_var = cases / pys^2
    )
  
  # collect the age-adjusted and age-specfic output into 1 vector
  out <- c(age_spec_ir_out, unlist(ir_out, use.names = TRUE))
  
  return(out)
}


# OLD versions - ignore! ----
# calculate_ir <- function(data, age_cat, fu_start, fu_end, event_flag, 
#                          std_pop, wt = NULL) {
#   
#   require(tidyverse)
#   # data <- temp
#   # age_cat <- c(59, 65, 70, 75, 80, 85, 120)
#   # fu_start <- "SURVEY_AGE"
#   # fu_end <- "MAIN_DEM_V1_END_AGE"
#   # event_flag <- "MAIN_DEM_V1_END_DEM_FLAG"
#   # wt <- "sw1"
#   
#   age_cat_labels <- c()
#   
#   for (i in 1:(length(age_cat) - 1)) {
#     # initialize the age interval
#     cat_start <- age_cat[i]
#     cat_end <- age_cat[i+1]
#     # create the variable names for py and case contribution
#     age_int <- paste0(cat_start, "_", cat_end)
#     vars <- paste0(c("py_", "case_"), age_int)
#     py_var <- vars[1]
#     case_var <- vars[2]
#     
#     # collect the age intervals for setting up the IR output table
#     age_cat_labels <- c(age_cat_labels, age_int)
#     
#     # calculate py and case contribution by age category
#     data <- data %>% 
#       mutate(
#         !! py_var := case_when(
#           # start of fu after age category or end of fu before age category: 
#           # contribute 0 pys
#           get(fu_start) > cat_end | get(fu_end) <= cat_start ~ 0, 
#           
#           # start of fu before age category and end of fu during age category: 
#           get(fu_start) <= cat_start & get(fu_end) <= cat_end ~ 
#             get(fu_end) - cat_start, 
#           
#           # start of fu before age category and end of fu after age category: 
#           get(fu_start) <= cat_start & get(fu_end) > cat_end ~ 
#             cat_end - cat_start, 
#           
#           # start and end of fu during age category:
#           get(fu_start) <= cat_end & get(fu_end) <= cat_end ~ 
#             get(fu_end) - get(fu_start), 
#           
#           # start of fu during age category and end of fu after age cateogory:
#           get(fu_start) <= cat_end & get(fu_end) > cat_end ~ 
#             cat_end - get(fu_start), 
#           
#           TRUE ~ NA_real_),  
#         
#         # if end of fu is during age category, case contribution is the dem flag. 
#         # otherwise, case contribution is zero
#         !! case_var := ifelse(
#           get(fu_end) <= cat_end & get(fu_end) > cat_start, 
#           get(event_flag), 0
#         )
#       )
#   }
#   
#   # using 2000 US census population as the standard population
#   age_spec_ir_out <- tibble(
#     age_range = age_cat_labels, 
#     pop_count = std_pop, 
#     pop_total = sum(pop_count), 
#     pop_prop = pop_count / pop_total, 
#     pys = NA, 
#     cases = NA, 
#     # crude_ir = NA,
#     # crude_ir_var = NA, 
#     # age_adj_ir = NA, 
#     # age_adj_ir_var = NA
#   )
#   
#   for (i in 1:length(age_spec_ir_out$age_range)) {
#     cat <- age_spec_ir_out$age_range[i]
#     age_spec_ir_out[i, "pys"] <- ifelse(
#       is.null(wt), sum(data[,paste0("py_", cat)]),
#       sum(data$sw1 * data[,paste0("py_", cat)])) 
#     age_spec_ir_out[i, "cases"] <- ifelse(
#       is.null(wt), sum(data[,paste0("case_", cat)]), 
#       sum(data$sw1 * data[,paste0("case_", cat)])
#     )
#   }
#   
#   age_spec_ir_out <- age_spec_ir_out %>% 
#     mutate(ir = cases / pys, 
#            ir_var = cases / pys^2, 
#            adj_ir = ir * pop_prop, 
#            adj_ir_var = ir_var * pop_prop^2) %>% 
#     select(-ir, -ir_var)
#   
#   ir_out <- age_spec_ir_out %>% 
#     summarise(across(pys:adj_ir_var, sum)) %>% 
#     mutate(
#       # generate crude ir and var
#       crude_ir = cases / pys,
#       crude_ir_var = cases / pys^2
#     )
#   
#   return(ir_out)
# }

# res <- calculate_ir(data = temp, age_cat = c(59, 65, 70, 75, 80, 85, 120), 
#                     fu_start = "SURVEY_AGE", fu_end = "MAIN_DEM_V1_END_AGE", 
#                     event_flag = "MAIN_DEM_V1_END_DEM_FLAG", 
#                     std_pop = pop, wt = "sw1")

