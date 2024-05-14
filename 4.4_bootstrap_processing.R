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

# analysis ----
# function to output point estimates and percentile CI's
boot_summary <- function(ethnicity) {
  # testing 
  # ethnicity <- "Chinese"
  
  # read in the bootstrap samples 
  boot <- list.files(
    paste0(path_to_box, path_to_proj, path_to_output, "bootstraps/"),
    pattern = ethnicity, full.names = TRUE
  ) %>% 
    lapply(read_csv) %>% 
    bind_rows() 
  
  # a helper function that converts IR into 1000 yrs and rounded to 2 digits
  py_round_fmt <- function(x) {
    formatC(round(x*1000, digits = 2), format = 'f', digits = 2)
  }

  # point estimate is the mean across 40 imputations using the full dataset (t = 0)
  point_est <- boot %>% 
    filter(t == 0) %>% # choose the full data results 
    select(-t, -imp, -warning, 
           -ends_with("cases"), -ends_with("pys"),
           -contains("AJ_time")) %>% 
    # IR by 1000 pys, rounded to 3 decimal places
    summarise_all(mean) %>% 
    mutate(across(-contains("AJ_est"), py_round_fmt))

  # collect the percentile CI's
  perc_CI <- boot %>%
    # pick bootstrap samples without error
    filter(warning == 0, t != 0) %>%
    # take the first 1000 bootstrap samples without error for each imputation
    group_by(imp) %>%
    slice_min(order_by = t, n = 1000) %>% 
    ungroup() %>% 
    select(-t, -imp, -warning, 
           -ends_with("cases"), -ends_with("pys"),
           -contains("AJ_time")) %>% 
    # find the 2.5 and 97.5 percentiles and convert to 1000 yrs
    summarise_all(quantile, probs = c(0.025, 0.975)) %>%
    mutate(across(-contains("AJ_est"), py_round_fmt))
  
  # prepare the plot data
  ## crude and age adjusted rates
  plot_data <- rbind(point_est, perc_CI) %>%
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
  
  # R1: cumulative risks
  plot_data_lifetime_risk <- rbind(point_est, perc_CI) %>%
    select(contains("AJ")) %>% 
    t() %>%
    as_tibble() %>% 
    rename(point = V1, lower = V2, upper = V3) %>%
    mutate(
      wt = rep(c("unweighted", "weighted"), each = 6),
      by_age = rep(seq(65, 90, by = 5), 2),
      across(c(point, lower, upper), function(x) round(x * 100, 0))
    )
  
  lifetime_risk <- plot_data_lifetime_risk %>% 
    mutate(est = paste0(point, " (", lower, ",", upper, ")")) %>% 
    pivot_wider(
      id_cols = c(by_age),
      names_from = wt, 
      values_from = est
    )

  return(list(table = out, 
              plot_data = plot_data, plot_data_age_spec = plot_data_age_spec,
              lifetime_risk = lifetime_risk,
              plot_data_lifetime_risk = plot_data_lifetime_risk
              ))
}

# run the function for each ethnicity
ch <- boot_summary("Chinese")
fi <- boot_summary("Filipino")
ja <- boot_summary("Japanese")
ko <- boot_summary("Korean")
pi <- boot_summary("Pacific Islander")
sa <- boot_summary("South Asian")
vi <- boot_summary("Vietnamese")
# mu <- boot_summary("Multiple Ethnicities")

# save the formatted tables

## incidence rates
results_table <- list(
  "Chinese" = ch$table,
  "Filipino" = fi$table,
  "Japanese" = ja$table,
  "Korean" = ko$table,
  "Pacific Islander" = pi$table,
  "South Asian" = sa$table,
  "Vietnamese" = vi$table
  # "Multiple Ethnicities" = mu$table
) %>% 
  bind_rows(.id = "ethnicity") %>% 
  pivot_longer(
    cols = c(unweighted, weighted), 
    names_to = "type", 
    values_to = "value"
  ) %>% 
  pivot_wider(
    id_cols = c(ethnicity, type),
    names_from = cat, 
    values_from = value
  ) %>% 
  arrange(type, ethnicity) 

## lifetime risk
lifetime_risk <- list(
  "Chinese" = ch$lifetime_risk,
  "Filipino" = fi$lifetime_risk,
  "Japanese" = ja$lifetime_risk,
  "Korean" = ko$lifetime_risk,
  "Pacific Islander" = pi$lifetime_risk,
  "South Asian" = sa$lifetime_risk,
  "Vietnamese" = vi$lifetime_risk
  # "Multiple Ethnicities" = mu$lifetime_risk
) %>% 
  bind_rows(.id = "ethnicity") %>% 
  pivot_longer(
    cols = c(unweighted, weighted), 
    names_to = "type", 
    values_to = "value"
  ) %>% 
  pivot_wider(
    id_cols = c(ethnicity, type),
    names_from = by_age, 
    values_from = value
  ) %>% 
  arrange(type, ethnicity) 

# write.xlsx( 
#   list(
#     "age_specific" = results_table %>% select(-crude, -`age adjusted`),
#     "age_adjusted" = results_table %>% 
#       select(ethnicity, type, crude, `age adjusted`) %>% 
#       pivot_wider(
#         id_cols = ethnicity,
#         names_from = type,
#         values_from = c(crude, `age adjusted`)
#       ),
#     "lifetime_risk" = lifetime_risk
#   ),
#   file = paste0(path_to_box, path_to_proj, path_to_output,
#                 "results_by_ethn.xlsx")
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
                     limits = c(0, 29)
                     ) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(
    # title = "Crude incidence rates by race/ethnicity and diabetes", 
    x = "Ethnicity", y = "Crude incidence rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )

# explore an option to include numerical labels on the figure
# plt_data %>% 
#   filter(type == "crude") %>% 
#   ggplot(aes(x = ethn, fill = ethn, 
#              alpha = wt, group = interaction(ethn, wt))) + 
#   geom_bar(aes(y = point), stat = "identity", position = position_dodge()) +
#   geom_errorbar(
#     stat = "identity",
#     aes(ymin = lower, ymax = upper), width = .2,
#     position = position_dodge(.9)
#   ) + 
#   # geom_text(
#   #   aes(x = ethn, y = point, label = round_pad(point, 1)),
#   #   alpha = 1, position = position_dodge(.9),
#   #   size = 3
#   # ) +
#   # geom_text(
#   #   aes(x = ethn, y = upper + 0.5, label = round_pad(upper, 1)),
#   #   alpha = 1, position = position_dodge(.9),
#   #   size = 3
#   # ) +
#   # geom_text(
# #   aes(x = ethn, y = lower - 0.5, label = round_pad(lower, 1)),
# #   alpha = 1, position = position_dodge(.9),
# #   size = 3
# # ) +
#   # geom_text(
#   #   aes(x = ethn, y = ifelse(wt == "unweighted", 27, 29), label = round_pad(point, 1)),
#   #   alpha = 1, position = position_dodge(.9),
#   #   size = 3
#   # ) +
#   # geom_text(
#   #   aes(x = ethn, y = ifelse(wt == "unweighted", 26, 28), label = paste0("(", round_pad(lower, 1), ",", round_pad(upper, 1), ")")),
#   #   alpha = 1, position = position_dodge(.9),
#   #   size = 3
#   # ) +
#   geom_text(
#     aes(x = ethn, y = upper + 1,
#         label = paste0(round_pad(point, 1), 
#         "\n(", round_pad(lower, 1), ",", round_pad(upper, 1), ")")),
#     alpha = 1, position = position_dodge(.9), 
#     size = 2.5
#   ) +
#   scale_y_continuous(breaks = seq(0, 27.5, by = 5), 
#                      limits = c(0, 29)
#   ) +
#   scale_alpha_discrete(range = c(0.5, 1)) +
#   labs(
#     # title = "Crude incidence rates by race/ethnicity and diabetes", 
#     x = "Ethnicity", y = "Crude incidence rates per 1,000 person-years") +
#   theme_bw() + 
#   theme(
#     # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#     legend.position = "none"
#   )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn_labelled3.png"),
#   width = 7, height = 5.5, units = "in"
# )

ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.png"),
  width = 7, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.tiff"),
  width = 7, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/crude_IR_x_ethn.pdf"),
  width = 7, height = 5.5, units = "in"
)


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
                     limits = c(0, 29)
  ) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(
    # title = "Age-adjusted incidence rates by race/ethnicity and diabetes", 
    x = "Ethnicity", y = "Age-standardized incidence rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "none"
  )

ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.tiff"),
  width = 7, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.png"),
  width = 7, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/age_adj_IR_x_ethn.pdf"),
  width = 7, height = 5.5, units = "in"
)

## age-specific rates ----

upper_limit <- 80
label_y <- 85

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
    
    # some upper limits are too large:
    # create a new upper limit
    upper_l = ifelse(upper < upper_limit, upper, NA),
    # specify y position for arrows 
    arrow_pos = ifelse(upper >= upper_limit, upper_limit, NA),
    # specify y position for labels
    label_pos = ifelse(upper >= upper_limit, label_y, NA)
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
    data = age_spec_plt_data %>% filter(less_than_5 == 0), 
    stat = "identity",
    aes(x = age_grp_num_new, ymin = lower, ymax = upper_l), width = 0.2,
  ) + 
  # add arrows for where the upper limits are outside the plot
  geom_segment(
    data = age_spec_plt_data %>% filter(less_than_5 == 0), 
    aes(x = age_grp_num_new, xend = age_grp_num_new, 
        y = lower, yend = arrow_pos),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = F
  ) + 
  # label the arrows with actual values of the upper limits
  geom_text(
    data = age_spec_plt_data %>% filter(less_than_5 == 0), 
    aes(x = age_grp_num_new, y = label_pos, 
        # label = paste0("UL=", round(upper, 2))
        label = round(upper, 0)),
    size = 2.5
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
  scale_y_continuous(limits = c(0, 88), breaks = seq(0, 80, by = 20)) + 
  labs(x = "Age group", y = "Age specific rates per 1,000 person-years") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )

# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_spec_IR_x_ethn.tiff"),
#   width = 10, height = 5.5, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_spec_IR_x_ethn.png"),
#   width = 10, height = 5.5, units = "in"
# )
# ggsave(
#   paste0(path_to_box, path_to_proj, "output/figures/age_spec_IR_x_ethn.pdf"),
#   width = 8, height = 5.5, units = "in"
# )


## lifetime risk at age 90 ----

lifetime_risk_plt_data <- list(
  "Chinese" = ch$plot_data_lifetime_risk,
  "Filipino" = fi$plot_data_lifetime_risk,
  "Japanese" = ja$plot_data_lifetime_risk,
  "Korean" = ko$plot_data_lifetime_risk,
  "Pacific Islander" = pi$plot_data_lifetime_risk,
  "South Asian" = sa$plot_data_lifetime_risk,
  "Vietnamese" = vi$plot_data_lifetime_risk
  # "Multiple Ethnicities" = mu$plot_data_lifetime_risk
) %>% 
  bind_rows(.id = "ethn") 

lifetime_risk_plt_data %>% 
  filter(by_age == 90) %>% 
  ggplot(aes(x = ethn, fill = ethn, 
             alpha = wt, group = interaction(ethn, wt))) + 
  geom_bar(aes(y = point), stat = "identity", position = position_dodge()) +
  geom_errorbar(
    stat = "identity",
    aes(ymin = lower, ymax = upper), width = .2,
    position = position_dodge(.9)
  ) + 
  scale_y_continuous(breaks = seq(0, 60, by = 10), 
                     limits = c(0, 60)
  ) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(
    x = "Ethnicity", y = "Cumulative risk (%)") +
  theme_bw() + 
  theme(
    legend.position = "none"
  )

ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/lifetime_risk_90_x_ethn.tiff"),
  width = 10, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/lifetime_risk_90_x_ethn.png"),
  width = 10, height = 5.5, units = "in"
)
ggsave(
  paste0(path_to_box, path_to_proj, "output/figures/lifetime_risk_90_x_ethn.pdf"),
  width = 7, height = 5.5, units = "in"
)


