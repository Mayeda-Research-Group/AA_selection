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
path_to_output <- "output/bootstrap_results/"


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

  return(list(table = out, plot_data = plot_data))
}

# run the function for each ethnicity
ch <- boot_summary("bootstrap_ir_Chinese_w_warn_index.csv")
fi <- boot_summary("bootstrap_ir_Filipino_w_warn_index.csv")
ja <- boot_summary("bootstrap_ir_Japanese_w_warn_index.csv")
ko <- boot_summary("bootstrap_ir_Korean_w_warn_index.csv")
pi <- boot_summary("bootstrap_ir_Pacific Islander_w_warn_index.csv")
sa <- boot_summary("bootstrap_ir_South Asian_w_warn_index_rerun_no_usborn.csv")
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
                     limits = c(0, 27.85)
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
#   width = 7, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.png"),
#   width = 7, units = "in"
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
                     limits = c(0, 27.85)
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
#   width = 7, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.png"),
#   width = 7, units = "in"
# )
