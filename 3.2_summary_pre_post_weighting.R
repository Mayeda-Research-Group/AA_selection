# generate figures and tables for manuscript: 
# (1) covariate balance plots by ethnicity before and after applying weight
# (2) summary statistics for propensity scores and weights

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "twang", "openxlsx", "reshape2", "survey",
       "gtsummary", "rcartocolor"
       # "mice", rlang"
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # , "lmtest", "ggpubr", "patchwork", "mitools", 
)

options(scipen = 999, digits = 8)

# set up paths and data ----
# path_to_box <- "C:/Users/Staff/Box/"
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
path_to_output <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/output/"

# load the saved objects containing covariate balance and weights
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/ch_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/fi_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/ja_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/ko_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/pi_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/sa_res.RData"))
load(file = paste0(path_to_box, path_to_datasets, "weighted_data/vi_res.RData"))


# new labels to be used for variables in covariance balance plot
labels <- c(
  "Age", "Female", "US born",
  "<= Grade 11", "Grade 12/high school diploma/GED", 
  "Technical/trade/vocational school/some college", "College degree",
  "Currently working", "Currently retired", 
  "Married or living as married", 
  "Living alone", "Household size = 2", "Household size = 3+", 
  "Household size-adjusted income", 
  "General health: excellent or very good", 
  "General health: good", 
  "General health: fair or poor", 
  "BMI", 
  "Ever had diabetes", 
  "Ever had hypertension", 
  "Smoking: never", "Smoking: former", "Smoking: current"
)

# order of their appearance in covariance balance plot;
# this matches with their order in Table 1
ordered_labels <- c(
  "Age", "Female", "US born",
  "<= Grade 11", "Grade 12/high school diploma/GED", 
  "Technical/trade/vocational school/some college", "College degree",
  "Married or living as married", 
  "Living alone", "Household size = 2", "Household size = 3+", 
  "Household size-adjusted income", 
  "Currently working", "Currently retired", 
  "General health: excellent or very good", 
  "General health: good", 
  "General health: fair or poor", 
  "BMI", 
  "Ever had diabetes", 
  "Ever had hypertension", 
  "Smoking: never", "Smoking: former", "Smoking: current"
)


# unweighted covariate balance ----

unw_cov_bal_data <- tibble(
  var = ch0_res$std_eff_sz$var, 
  "Chinese" = ch0_res$std_eff_sz$mean, 
  "Filipino" = fi0_res$std_eff_sz$mean, 
  "Japanese" = ja0_res$std_eff_sz$mean, 
  "Korean" = ko0_res$std_eff_sz$mean, 
  "Pacific Islander" = pi0_res$std_eff_sz$mean, 
  "South Asian" = sa0_res$std_eff_sz$mean, 
  "Vietnamese" = vi0_res$std_eff_sz$mean
  ) %>% 
  filter(
    # remove some redundant rows for simplicity
    str_detect(var, "No|Male", negate = TRUE)
  ) %>% 
  mutate(
    # label = labels, 
    new_label = factor(labels, levels = ordered_labels))


# weighted covariate balance ----

wt_cov_bal_data <- tibble(
  var = ch3_res$std_eff_sz$var, 
  "Chinese" = ch3_res$std_eff_sz$mean, 
  "Filipino" = fi5_res$std_eff_sz$mean, 
  "Japanese" = ja3_res$std_eff_sz$mean, 
  "Korean" = ko3.2_res$std_eff_sz$mean, 
  "Pacific Islander" = pi2_res$std_eff_sz$mean, 
  "South Asian" = sa5_res$std_eff_sz$mean, 
  "Vietnamese" = vi4_res$std_eff_sz$mean
) %>% 
  filter(
    # remove some redundant rows for simplicity
    str_detect(var, "No|Male", negate = TRUE)
  ) %>% 
  mutate(
    # label = labels, 
    new_label = factor(labels, levels = ordered_labels))


# combine the unweighted and weighted plots as two layers
cov_bal_plot <- ggplot() +
  geom_point(data = melt(unw_cov_bal_data), 
             aes(x = new_label, y = value, color = variable), 
             shape = "circle open") + 
  geom_point(data = melt(wt_cov_bal_data), 
             aes(x = new_label, y = value, color = variable), 
             shape = "circle") + 
  geom_hline(yintercept = 0) +
  ylim(-1.25, 1.25) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  labs(
    y = expression(
      atop(
        "Standardized mean difference",
        paste("[mean(", KPNC[k], ") - mean(", CHIS[k], ")] / SD(", CHIS[k], ")")
      )
    ), 
    x = NULL, 
    color = "Race/Ethnicity"
  ) +
  theme_bw() +
  theme(
    # axis.text.y = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_line(linetype = c("solid", "dashed")),
    # aspect.ratio = 5 / 3
    aspect.ratio = 9 / 5
  ) + 
  facet_wrap( ~ variable, ncol = 4) +
  theme(legend.position = "none")

ggsave(
  plot = cov_bal_plot,
  file = paste0(path_to_box, path_to_output, "figures/cov_bal_unw.png"),
  width = 8, units = "in"
)

# rotate the figure for submission (landscape)
ggsave(
  plot = cov_bal_plot,
  file = paste0(path_to_box, path_to_output, "figures/cov_bal_unw.pdf"),
  width = 8, height = 6, units = "in"
)

# summary after weighting ----

## numerical summaries ----

# function to obtain summary statistics for propensity score or weights
# obtained from using CHIS and the last copy of imputed RPGEH
ps_wt_summarize <- function(data, ethn_str, var) {
  # testing
  # data <- ch3_res$wt_df
  # ethn_str <- "Chinese"
  # var <- "p1"
  
  if (var == "p1") {
    # for propensity scores, we want weighted summary statistics for CHIS
    # recall that smplwt is 1 for RPGEH and survey weights for CHIS
    dsn <- data %>% 
      # obtain CHIS and the last copy of imputed RPGEH
      filter(imp %in% c(0, 40)) %>% 
      svydesign(~ 1, data = ., weights = ~ smplwt)
    out <- tbl_svysummary(
      data = dsn, 
      by = in_chis, 
      statistic = all_continuous() ~ c("{mean}", "{sd}", "{min}", "{p10}", 
                                       "{median}", "{p90}", "{max}"), 
      include = c(p1),
      type = list(p1 ~ "continuous2"),
      digits = all_continuous() ~ 4
    ) %>% as_tibble()
    
    out <- rbind(
      list("Ethnicity", ethn_str, ethn_str), 
      list("Sample", "Kaiser", "CHIS"), 
      out[-1,]
    )
  
  } else {
    # for weights, we only want summary statistics for RPGEH
    out <- data %>% 
      filter(imp == 40) %>% # choose the last copy of RPGEH
      summarise(mean = mean(sw1), 
                sd = sd(sw1),
                min = min(sw1), 
                quantile_10 = quantile(sw1, probs = 0.1), 
                median = median(sw1), 
                quantile_90 = quantile(sw1, probs = 0.9), 
                max = max(sw1)) %>% 
      round(digits = 3)
    out <- rbind(ethn_str, t(out))
  }
  
  return(out)
}


# collect the summary statistics and save as excel
ps_stats <- cbind(
  ps_wt_summarize(ch3_res$wt_df, "Chinese", "p1"), 
  ps_wt_summarize(fi5_res$wt_df, "Filipino", "p1")[, -1], 
  ps_wt_summarize(ja3_res$wt_df, "Japanese", "p1")[, -1], 
  ps_wt_summarize(ko3.2_res$wt_df, "Korean", "p1")[, -1], 
  ps_wt_summarize(pi2_res$wt_df, "Pacific Islander", "p1")[, -1], 
  ps_wt_summarize(sa5_res$wt_df, "South Asian", "p1")[, -1], 
  ps_wt_summarize(vi4_res$wt_df, "Vietnamese", "p1")[, -1]
)

wt_stats <- cbind(
  ps_wt_summarize(ch3_res$wt_df, "Chinese", "sw1"), 
  ps_wt_summarize(fi5_res$wt_df, "Filipino", "sw1"), 
  ps_wt_summarize(ja3_res$wt_df, "Japanese", "sw1"), 
  ps_wt_summarize(ko3.2_res$wt_df, "Korean", "sw1"), 
  ps_wt_summarize(pi2_res$wt_df, "Pacific Islander", "sw1"), 
  ps_wt_summarize(sa5_res$wt_df, "South Asian", "sw1"), 
  ps_wt_summarize(vi4_res$wt_df, "Vietnamese", "sw1")
)

wt_stats <- cbind(rownames(wt_stats), wt_stats)

# write.xlsx(list(ps = ps_stats, wt = wt_stats),
#            file = paste0(path_to_box, path_to_output,
#                          "propensity_score_weights_summary.xlsx"))

## plots ----
# collect propensity scores and weights for CHIS and the last copy of RPGEH
ps_wt <- rbind(
  ch3_res$wt_df %>% filter(imp %in% c(0, 40)), 
  fi5_res$wt_df %>% filter(imp %in% c(0, 40)), 
  ja3_res$wt_df %>% filter(imp %in% c(0, 40)), 
  ko3.2_res$wt_df %>% filter(imp %in% c(0, 40)), 
  pi2_res$wt_df %>% filter(imp %in% c(0, 40)), 
  sa5_res$wt_df %>% filter(imp %in% c(0, 40)), 
  vi4_res$wt_df %>% filter(imp %in% c(0, 40))
) %>% 
  mutate(study = ifelse(in_chis == 1, "CHIS", "Kaiser"))

### propensity score ----
# weighted by sampling weights
ps_plot <- ggplot(ps_wt, aes(x = p1, group = study, color = study)) +
  # note that we apply survey weights here for CHIS and weight = 1 for RPGEH
  geom_density(aes(fill = study, weight = smplwt), alpha = 0.2) +
  xlim(0, 0.3) +
  labs(y = "Density", 
       x = "Propensity Score (axis truncated at 0.3 for readability)", 
       color = "Study") +
  scale_color_grey(start = 0, end = 0.5) +
  scale_fill_grey(start = 0, end = 0.5) +
  # scale_color_carto_d(palette = "Vivid", direction = -1) +
  # scale_fill_carto_d(palette = "Vivid", direction = -1) + 
  guides(fill = "none") +
  theme_bw() +
  facet_wrap(~ H_ethn, ncol = 2)

# ggsave(
#   plot = ps_plot,
#   file = paste0(path_to_box, path_to_output, "figures/propensity_score.png"),
#   width = 7, height = 6.5, units = "in"
# )
  

### weights ----
wt_plot <- ps_wt %>% filter(study == "Kaiser") %>% 
  ggplot(aes(x = sw1, group = H_ethn, color = study)) +
  geom_density(aes(fill = study), alpha = 0.2) +
  xlim(0, 10) +
  labs(y = "Density", 
       x = "Stabilized inverse odds of selection weights\n(axis truncated at 10 for readability)") +
  # scale_color_carto_d(palette = "Vivid") +
  scale_color_grey(start = 0.5) +
  scale_fill_grey(start = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ H_ethn, ncol = 2)

# ggsave(
#   plot = wt_plot,
#   file = paste0(path_to_box, path_to_output, "figures/weights.png"),
#   width = 7, height = 6.5, units = "in"
# )


  