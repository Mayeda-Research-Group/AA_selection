# process bootstrap results

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "openxlsx"
       # "mice", rlang", "twang", 
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "boot", "lmtest", "mitools", "ggpubr", "patchwork", "RColorBrewer"
)

options(scipen = 999, digits = 8)

# set up paths and data ----

path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_proj <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/"
path_to_output <- "output/bootstrap_results_strat/" 
# modified on Mar 23: process bootstrap results using stratified resampling

# analysis ----
# function to output point estimates and percentile CI's
boot_summary <- function(data_path) {
  # testing 
  # data_path <- "bootstrap_ir_Chinese_w_warn_index.csv"
  
  # read in the bootstrap samples by data path
  boot <- read_csv(
    paste0(path_to_box, path_to_proj, path_to_output, 
           "bootstrap_samples/", data_path)
  )
  
  # a helper function that converts IR into 1000 yrs and rounded to 2 digits
  py_round_fmt <- function(x) {
    formatC(round(x*1000, digits = 2), format = 'f', digits = 2)
  }
  
  # NEW from 01/26/2023: fix age-specific rates
  # The age specific rates shown in these bootstrap samples are after being 
  # multiplied with the proportion of the corresponding age category in the
  # standard population to calculate age-adjusted IRs, but to report the actual
  # age specific rates, this multiplication step is not needed. Here, I am 
  # reversing this step. 

  std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
  std_pop_prop <- std_pop / sum(std_pop)
  boot <- boot %>% 
    mutate(
      across(ends_with("59_65"), function(x) {x / std_pop_prop[1]}),
      across(ends_with("65_70"), function(x) {x / std_pop_prop[2]}),
      across(ends_with("70_75"), function(x) {x / std_pop_prop[3]}),
      across(ends_with("75_80"), function(x) {x / std_pop_prop[4]}),
      across(ends_with("80_85"), function(x) {x / std_pop_prop[5]}),
      across(ends_with("85_120"), function(x) {x / std_pop_prop[6]}),
    )
  
  
  # point estimate is the mean across 40 imputations using the full dataset (t = 0)
  point_est <- boot %>% 
    filter(t == 0) %>% # choose the full data results 
    select(-t, -imp, -warning, -ends_with("cases"), -ends_with("pys")) %>% 
    # IR by 1000 pys, rounded to 3 decimal places
    summarise_all(mean) %>% 
    mutate_all(py_round_fmt) %>% 
    as.matrix()

  # collect the percentile CI's
  perc_CI <- boot %>%
    # pick bootstrap samples without error
    filter(warning == 0, t != 0) %>%
    # take the first 1000 bootstrap samples without error for each imputation
    group_by(imp) %>%
    slice_min(order_by = t, n = 1000) %>% 
    ungroup() %>% 
    select(-t, -imp, -warning, -ends_with("cases"), -ends_with("pys")) %>% 
    # find the 2.5 and 97.5 percentiles and convert to 1000 yrs
    summarise_all(quantile, probs = c(0.025, 0.975)) %>%
    mutate_all(py_round_fmt) %>% 
    as.matrix()
  
  # prepare the plot data
  ## crude and age adjusted rates
  plot_data <- rbind(point_est, perc_CI) %>%
    as_tibble() %>% 
    select(ends_with("ir")) %>% 
    t() %>%
    as_tibble() %>%
    rename(point = V1, lower = V2, upper = V3) %>%
    mutate(
      wt = rep(c("unweighted", "weighted"), 2), 
      type = rep(c("age adjusted", "crude"), each = 2)
    ) 
  
  ## age specific rates
  plot_data_age_spec <- rbind(point_est, perc_CI) %>%
    as_tibble() %>% 
    select(contains("spec")) %>% 
    t() %>%
    as_tibble() %>% 
    rename(point = V1, lower = V2, upper = V3) %>%
    mutate(
      wt = rep(c("unweighted", "weighted"), each = 6),
      age_grp = rep(c("60-65", "65-70", "70-75", "75-80", "80-85", "85+"), 2)
    )

  # prepare the output table
  res <- paste0(
    point_est[1, ], 
    " (", perc_CI[1, ], ",", perc_CI[2, ], ")"
  )
  
  out <- tibble(
    cat = c("60-65", "65-70", "70-75", "75-80", "80-85", "85+", 
            "crude", "age adjusted"), 
    unweighted = c(res[1:6], res[c(15,13)]), 
    weighted = c(res[7:12], res[c(16,14)]), 
  )

  return(list(table = out, plot_data = plot_data, plot_data_age_spec = plot_data_age_spec))
}

# run the function for each ethnicity
ch <- boot_summary("bootstrap_ir_Chinese_w_warn_index.csv")
fi <- boot_summary("bootstrap_ir_Filipino_w_warn_index.csv")
ja <- boot_summary("bootstrap_ir_Japanese_w_warn_index.csv")
ko <- boot_summary("bootstrap_ir_Korean_w_warn_index.csv")
pi <- boot_summary("bootstrap_ir_Pacific Islander_w_warn_index.csv")
sa <- boot_summary("bootstrap_ir_South Asian_w_warn_index_SA_rerun.csv")
vi <- boot_summary("bootstrap_ir_Vietnamese_w_warn_index.csv")
mu <- boot_summary("bootstrap_ir_Multiple Ethnicities_w_warn_index.csv")

# save the formatted table as output
# write.xlsx(
#   list(
#     "Chinese" = ch$table,
#     "Filipino" = fi$table,
#     "Japanese" = ja$table,
#     "Korean" = ko$table,
#     "Pacific Islander" = pi$table,
#     "South Asian" = sa$table,
#     "Vietnamese" = vi$table,
#     "Multiple Ethnicities" = mu$table
#   ),
#   file = paste0(path_to_box, path_to_proj, path_to_output,
#                 "IR_results_by_ethn.xlsx")
# )

# figures ----

## crude and age-adjusted rates ----

plt_data <- rbind(
  ch$plot_data, fi$plot_data, ja$plot_data, ko$plot_data, 
  pi$plot_data, sa$plot_data, vi$plot_data) %>% 
  mutate(
    ethn = rep(c("Chinese", "Filipino", "Japanese", "Korean", "Pacific Islander", 
                 "South Asian", "Vietnamese"), each = 4)) %>% 
  mutate(
    point = as.numeric(point), 
    lower = as.numeric(lower), 
    upper = as.numeric(upper), 
    # # some upper limits are too large; replace with a limit for ggplot
    # upper_l = ifelse(upper >= 25, 25, upper)
  )


plt_data %>% 
  filter(type == "crude") %>% 
  ggplot(aes(x = ethn, fill = ethn, 
             alpha = wt, group = interaction(ethn, wt))) + 
  geom_bar(aes(y = point), stat = "identity", position = position_dodge()) +
  geom_errorbar(
    stat = "identity",
    aes(ymin = lower, ymax = upper), width = .2,
    position = position_dodge(.9)
    ) + 
  scale_y_continuous(breaks = seq(0, 27.5, by = 5), 
                     limits = c(0, 28)
                     ) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(
    # title = "Crude incidence rates by race/ethnicity and diabetes", 
    x = "Race/Ethnicity", y = "Crude incidence rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )

# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.tiff"),
#   width = 7, height = 5.5, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.png"),
#   width = 7, height = 5.5, units = "in"
# )

plt_data %>% 
  filter(type == "age adjusted") %>% 
  ggplot(aes(x = ethn, fill = ethn, 
             alpha = wt, group = interaction(ethn, wt))) + 
  geom_bar(aes(y = point), stat = "identity", position = position_dodge()) +
  geom_errorbar(
    stat = "identity",
    aes(ymin = lower, ymax = upper), width = .2,
    position = position_dodge(.9)
  ) + 
  scale_y_continuous(breaks = seq(0, 27.5, by = 5), 
                     limits = c(0, 28)
  ) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(
    # title = "Age-adjusted incidence rates by race/ethnicity and diabetes", 
    x = "Race/Ethnicity", y = "Age-adjusted incidence rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "none"
  )

# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.tiff"),
#   width = 7, height = 5.5, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.png"),
#   width = 7, height = 5.5, units = "in"
# )

## age-specific rates ----

upper_limit <- 80

age_spec_plt_data <- rbind(
  ch$plot_data_age_spec, fi$plot_data_age_spec, ja$plot_data_age_spec, 
  ko$plot_data_age_spec, pi$plot_data_age_spec, sa$plot_data_age_spec, 
  vi$plot_data_age_spec
) %>% 
  mutate(
    # add in ethnicity and age group by hand
    ethn = rep(c("Chinese", "Filipino", "Japanese", "Korean", "Pacific Islander", 
                 "South Asian", "Vietnamese"), each = 12),
    age_grp_num = rep(1:6, 14), 
    
    # stagger the age group to plot unw and wt side by side
    age_grp_num_new = ifelse(wt == "unweighted", age_grp_num - 0.2, age_grp_num + 0.2), 
    
    # add indicator for whether there were less than 5 events for each cell
    less_than_5 = case_when(
      age_grp == "60-65" & ethn != "Filipino"  ~ 1, 
      age_grp == "65-70" & ethn %in% c("Korean", "Pacific Islander", 
                                       "South Asian", "Vietnamese") ~ 1, 
      age_grp == "70-75" & ethn == "Vietnamese" ~ 1, 
      TRUE ~ 0
    ),
    
    # format the numbers
    point = as.numeric(point), 
    lower = as.numeric(lower), 
    upper = as.numeric(upper), 
    # some upper limits are too large: create an indicator for 
    # whether it exceeds a given value
    # upper_l = ifelse(upper >= 80, 80, upper), 
    upper_ind = ifelse(upper >= upper_limit, 1, 0)
  )


ggplot() + 
  # for where there are more than 5 events, plot the bars
  geom_bar(
    data = age_spec_plt_data %>% filter(less_than_5 == 0), 
    aes(x = age_grp_num_new, y = point, fill = ethn, alpha = wt),
    stat = "identity", width = 0.4
  ) + 
  # add error bars for where the upper limits are contained within the plot
  geom_errorbar(
    data = age_spec_plt_data %>% filter(upper_ind == 0, less_than_5 == 0), 
    stat = "identity",
    aes(x = age_grp_num_new, ymin = lower, ymax = upper), width = .2,
  ) + 
  # for where the upper limits are not contained within the plot, 
  # first add in the bars for lower limits
  geom_errorbar(
    data = age_spec_plt_data %>% filter(upper_ind == 1),
    stat = "identity",
    aes(x = age_grp_num_new, ymax = upper_limit + 5, ymin = lower), width = .2,
  ) + 
  # then draw the arrows from the lower limits to the top of the plot
  geom_segment(
    data = age_spec_plt_data %>% filter(upper_ind == 1),
    aes(x = age_grp_num_new, xend = age_grp_num_new, 
        y = lower, yend = upper_limit), 
    arrow = arrow(length = unit(0.2, "cm"))
  ) + 
  # for where there are less than 5 events, put a cross
  geom_point(
    data = age_spec_plt_data %>% filter(less_than_5 == 1, wt == "unweighted"), 
    aes(x = age_grp_num, y = 0), shape = 4
  ) + 
  facet_wrap(~ ethn) + 
  scale_alpha_discrete(range = c(0.5, 1)) + 
  scale_x_continuous(
    breaks = 1:6, minor_breaks = NULL, 
    labels = c("60-65", "65-70", "70-75", "75-80", "80-85", "85+")
  ) + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20)) + 
  labs(x = "Age group", y = "Age specific rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )

ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/age_spec_IR_x_ethn.tiff"),
  width = 7, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/age_spec_IR_x_ethn.png"),
  width = 7, height = 5.5, units = "in"
)


