# this script contains sensitivity analysis on how to include age in the 
# GLM weighting model 
# we explore 2 options: linear and cubic bspline
# we examine covariate balance by std mean diff and weighted density plots, and 
# finally calculate and compare IR's using the two sets of weights 
# this sensitivity analysis is done on one imputation of the RPGEH dataset
# as a demonstration 

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "twang", "ggpubr", "splines", "rlang"
       # "mice", rlang", "RColorBrewer", "mitools", "patchwork", 
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "openxlsx", "lmtest"
)

options(scipen = 999, digits = 8)

# set up paths and data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
path_to_output <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/output/"

# load the prepped data generated in 4.1 and the IR calculation function
load(paste0(path_to_box, path_to_datasets, "bootstrap_prep_data.RData"))
source(paste0(path_to_box,
              "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
              "cleaned_scripts/function_IR_calc.R"))

# prepare datasets and other objects ----
# subset to one copy of the imputation
chis_rpgeh_imp1 <- chis_rpgeh %>% filter(imp %in% c(0, 1)) 
# %>% mutate(H_age_cat = cut(H_age, breaks = c(59, 65, 70, 75, 80, 85, 90)))

# calculate unconditional RPGEH sampling odds 
# overall 
p_rpgeh <- sum(chis_rpgeh_imp1$in_rpgeh) / sum(chis_rpgeh_imp1$smplwt)
uncon_rpgeh_odds <- p_rpgeh / (1 - p_rpgeh)
# by ethnicity
p_rpgeh_ethn <- chis_rpgeh_imp1 %>% group_by(H_ethn) %>%
  summarize(marginal_p = sum(in_rpgeh) / sum(smplwt),
            odds = marginal_p / (1 - marginal_p)) 


# append to the glm formulas to include a bspline age term
fmls_bspline <- lapply(
  fmls, 
  function(x){
    list("linear + other vars" = x, 
         "bspline + other vars" = str_replace(x, "H_age", "bs(H_age)"))
  })

# vars to check balance for
harmonized_norace <- c(
  "H_age", "H_female", "H_usborn", "H_edu_4", "H_work", "H_retired", 
  "H_marit", "H_hhsize_3", "H_income_pp", "H_health_3", "H_bmi", 
  "H_diab", "H_hyp", "H_smk"
)

# nicer labels 
labels <- c(
  "Age", "Female", "US born",
  "<= Grade 11", "Grade 12/HS diploma/GED", 
  "Technical/trade/vocational school/some college", "College and above",
  "Currently working", "Currently retired", 
  "Married or living as married", 
  "Living alone", "Household size = 2", "Household size = 3+", 
  "Household-adjusted income", 
  "General health: excellent or very good", 
  "General health: good", 
  "General health: fair or poor", 
  "BMI", 
  "Ever had diabetes", 
  "Ever had hypertension", 
  "Smoking: never", "Smoking: former", "Smoking: current"
)

ordered_labels <- c(
  "Age", "Female", "US born",
  "<= Grade 11", "Grade 12/HS diploma/GED", 
  "Technical/trade/vocational school/some college", "College and above",
  "Married or living as married", 
  "Living alone", "Household size = 2", "Household size = 3+", 
  "Household-adjusted income", 
  "Currently working", "Currently retired", 
  "General health: excellent or very good", 
  "General health: good", 
  "General health: fair or poor", 
  "BMI", 
  "Ever had diabetes", 
  "Ever had hypertension", 
  "Smoking: never", "Smoking: former", "Smoking: current"
)

# check covariate balance and generate weights given fmls ----
# the following function is a similar version to what's in 3.1
# only much simplified

wt_cov_bal <- function(data, ethn, fmls = NULL) {
  # testing
  # data <- chis_rpgeh_imp1
  # ethn <- "Chinese"
  # fmls <- fmls_bspline$Chinese # a list of formulas to iterate over
  
  # subset to a specific ethnicity
  ethn_subset <- data %>% filter(H_ethn == ethn)
  
  # check unweighted cov bal
  temp <- ethn_subset %>%
    as.data.frame() %>%
    bal.stat(
      data = .,
      vars = harmonized_norace,
      w.all = .[, "smplwt"],
      treat.var = 'in_chis',
      sampw = .$smplwt,
      estimand = "ATT",
      multinom = F
    )
  
  std_eff_sz <- tibble(
    var = rownames(temp$results),
    unw = -temp$results$std.eff.sz # reverse the sign
  )
  
  # weighted: iterate over the list of fmls for the weighting model
  for (i in 1:length(fmls)) {
    # i <- 1
    # fit glm model given each formula
    fml <- fmls[i]
    mod <- ethn_subset %>% 
      glm(
        as.formula(unlist(fml)),
        data = .,
        family = binomial(link = "logit"),
        weights = .$smplwt
      )
    
    # calculate propensity scores and weights
    ethn_subset[, paste0("p", i)] <- predict.glm(mod, type = "response")
    rpgeh_odds <- p_rpgeh_ethn %>% filter(H_ethn == ethn) %>% pull(odds)
    ethn_subset[, paste0("iow", i)] <- (1 - ethn_subset[, paste0("p", i)]) / 
      ethn_subset[, paste0("p", i)]
    ethn_subset[, paste0("sw", i)] <- 1
    ethn_subset[ethn_subset$in_rpgeh == 1, paste0("sw", i)] <- rpgeh_odds * 
      ethn_subset[ethn_subset$in_rpgeh == 1, paste0("iow", i)]
    ethn_subset[, paste0("sw", i, "_twang")] <- ethn_subset[, paste0("sw", i)] *
      ethn_subset[, "smplwt"] 
    # equivalent to the following, but I have more than 1 set of weights here
    # ethn_subset %>% mutate(
    #   iow1 = (1 - p1) / p1,
    #   sw1 = ifelse(in_rpgeh == 1, iow1 * rpgeh_odds, 1),
    #   sw1_twang = sw1 * smplwt
    # ) %>% View()
    
    # check covariate balance after weighting
    temp <- ethn_subset %>%
      as.data.frame() %>%
      bal.stat(
        data = .,
        vars = harmonized_norace,
        w.all = .[, paste0("sw", i, "_twang")],
        treat.var = 'in_chis',
        sampw = .$smplwt,
        estimand = "ATT",
        multinom = F
      )
    
    # append the std effect sizes
    std_eff_sz[, paste0("wt", i)] <- -temp$results$std.eff.sz # reverse the sign
    
  }
  
  # plot age densities in CHIS and in RPGEH before and after weighting
  plot_list <- vector(mode = "list", length = length(fmls))
  for (i in 1:length(fmls)) {
    plot_list[[i]] <-
      ggplot() +
      geom_density(
        data = filter(chis_rpgeh_imp1, in_chis == 1, H_ethn == ethn),
        aes(x = H_age, linetype = "CHIS: survey-weighted age", color = "CHIS", 
            weight = smplwt), 
      ) + 
      geom_density(
        data = filter(ethn_subset, in_chis == 0),
        aes(x = H_age, 
            linetype = "RPGEH: harmonized age top-coded at 85", 
            color = "unweighted RPGEH"
        )
      ) +
      geom_density(
        data = filter(ethn_subset, in_chis == 0),
        aes(
          x = H_age,
          linetype = "RPGEH: weighted harmonized age",
          weight = .data[[paste0("sw", i)]],
          color = "weighted RPGEH"
        )
      ) + 
      scale_linetype_manual(
        "Age variables",
        values = c("CHIS: survey-weighted age" = "solid",
                   "RPGEH: harmonized age top-coded at 85" = "dotted", 
                   "RPGEH: weighted harmonized age" = "dashed")
      ) + 
      scale_color_manual(
        "Study", 
        values = c("CHIS" = "black", "unweighted RPGEH" = "grey70",
                   "weighted RPGEH" = "grey50"), 
        guide = "none"
      ) + 
      labs(x = "age") + 
      scale_x_continuous(breaks = seq(60, 90, 10)
                         # , limits = c(60, 85)
      ) +
      labs(subtitle = paste0(ethn, ": weighted with age - ", names(fmls)[i])) +
      theme_bw() 
  }
  # arrange the plots
  age.dist <- ggarrange(plotlist = plot_list,
                        nrow = length(fmls), ncol = 1,
                        common.legend = TRUE,
                        legend = "right"
  )
  
  return(list(
    wt_df = ethn_subset,
    std_eff_sz = std_eff_sz,
    age_dist_plot = age.dist
  ))
  
}

# run the wt_cov_bal function ----
ch_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Chinese", 
                           fmls = fmls_bspline$Chinese)
fi_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Filipino", 
                           fmls = fmls_bspline$Filipino)
ja_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Japanese", 
                           fmls = fmls_bspline$Japanese)
ko_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Korean", 
                           fmls = fmls_bspline$Korean)
pi_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Pacific Islander", 
                           fmls = fmls_bspline$`Pacific Islander`)
sa_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "South Asian", 
                           fmls = fmls_bspline$`South Asian`)
vi_lin_vs_bs <- wt_cov_bal(data = chis_rpgeh_imp1, ethn = "Vietnamese", 
                           fmls = fmls_bspline$Vietnamese)

# save the density plots
ethn_abbr <- c("ch", "fi", "ja", "ko", "pi", "sa", "vi")

for (ethn_item in ethn_abbr) {
  ggsave(
    filename = paste0(
      path_to_box, path_to_output, "figures/linear vs bspline age/", 
      "age_density_", ethn_item, ".png"),
    plot = eval(parse_expr(paste0(ethn_item, "_lin_vs_bs$age_dist_plot"))),
    bg = "white", height = 7, units = "in")
}

# collect the std mean diff and plot by ethnicity
overall_std_eff_sz <- rbind(
  cbind(ch_lin_vs_bs$std_eff_sz, ethn = "Chinese"), 
  cbind(fi_lin_vs_bs$std_eff_sz, ethn = "Filipino"), 
  cbind(ja_lin_vs_bs$std_eff_sz, ethn = "Japanese"), 
  cbind(ko_lin_vs_bs$std_eff_sz, ethn = "Korean"), 
  cbind(pi_lin_vs_bs$std_eff_sz, ethn = "Pacific Islander"), 
  cbind(sa_lin_vs_bs$std_eff_sz, ethn = "South Asian"), 
  cbind(vi_lin_vs_bs$std_eff_sz, ethn = "Vietnamese")
)

overall_std_eff_sz %>% 
  filter(
    # remove some redundant rows for simplicity
    str_detect(var, "No|Male", negate = TRUE)
  ) %>%
  mutate(new_label = factor(rep(labels, 7), levels = ordered_labels)) %>% 
  ggplot(., aes(x = new_label)) +
  # geom_point(aes(y = unw, shape = "unweighted")) + 
  geom_point(aes(y = wt1, shape = "weighted with linear age")) + 
  geom_point(aes(y = wt2, shape = "weighted with bspline age")) + 
  scale_shape_manual(
    NULL,
    values = c(
      # "unweighted" = 20, 
      "weighted with linear age" = 20,
      "weighted with bspline age" = 2)
  ) + 
  geom_hline(yintercept = 0) +
  labs(
    y = "Standardized mean difference\n(RPGEH-CHIS)/SD[CHIS]",
    x = NULL) +
  theme_bw() +
  theme(
    # axis.title.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    # panel.grid.major.y = element_line(colour = c("grey85", "grey95")),
    panel.grid.major.y = element_line(linetype = c("solid", "dashed")),
    # legend.text = element_text(size = 7), 
    # aspect.ratio = 5 / 3
  ) + 
  ylim(-1, 1.25) +
  scale_x_discrete(limits = rev) +
  coord_flip() + 
  facet_wrap( ~ ethn)

# save the plot
ggsave(
  filename = paste0(
    path_to_box, path_to_output, "figures/linear vs bspline age/", 
    "cov_bal_linear_vs_bspline_age.png"),
  plot = last_plot(),
  height = 8, units = "in"
)

# calculate IR ----
# a function that calculates the IR's given the two types of weights 
compare_ir <- function(data, ethn) {
  # data <- ch_lin_vs_bs$wt_df
  # ethn <- "Chinese"
  boot_data <- data %>% 
    filter(in_rpgeh == 1) %>% # only need RPGEH participants
    select(ID, H_ethn, sw1, sw2) %>% # only need ethn and weights
    left_join(., pys_data, "ID") # join in tte data
  
  lin_ir <- boot_data %>%
    calculate_ir(
      data = ., 
      age_cat = age_cat_labels, 
      # this is the population we standardize to 
      std_pop = c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
    ) %>% 
    enframe() %>% 
    mutate(H_ethn = ethn, 
           age_term = "linear")
  
  bs_ir <- boot_data %>%
    # because of how calculate_ir() is set up, 
    # I have to rename the bspline weights
    select(-sw1) %>% 
    rename(sw1 = sw2) %>% 
    calculate_ir(
      data = ., 
      age_cat = age_cat_labels, 
      # this is the population we standardize to 
      std_pop = c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
    ) %>% 
    enframe() %>% 
    mutate(H_ethn = ethn, 
           age_term = "bspline")
  
  return(rbind(lin_ir, bs_ir))
}


ir_by_ethn <- rbind(
  compare_ir(ch_lin_vs_bs$wt_df, "Chinese"), 
  compare_ir(fi_lin_vs_bs$wt_df, "Filipino"),
  compare_ir(ja_lin_vs_bs$wt_df, "Japanese"),
  compare_ir(ko_lin_vs_bs$wt_df, "Korean"),
  compare_ir(pi_lin_vs_bs$wt_df, "Pacific Islander"),
  compare_ir(sa_lin_vs_bs$wt_df, "South Asian"),
  compare_ir(vi_lin_vs_bs$wt_df, "Vietnamese")
) %>% 
  # only want to look at aggregate IR's 
  filter(name %in% c("unw_adj_ir", "wt_adj_ir", "unw_crude_ir", "wt_crude_ir")) %>% 
  # create some new variables for plot purposes
  mutate(
    value = 1000*value, # IR per 1000 person-years
    wt = ifelse(str_detect(name, "unw"), "unweighted", "weighted"),
    age_adj = ifelse(str_detect(name, "adj"), "age-adjusted", "crude") %>% 
      factor(levels = c("crude", "age-adjusted")), 
    wt_age_term = ifelse(wt == "weighted", age_term, "unweighted") %>% 
      factor(levels = c("unweighted", "linear", "bspline"))
  )

ir_by_ethn %>% 
  # unweighted results are the same regardless of age terms
  # so I am removing one of them for redundancy
  filter(!(wt == "unweighted" & age_term == "bspline")) %>% 
  ggplot(., aes(x = value)) +
  geom_bar(aes(y = wt_age_term, fill = H_ethn, alpha = wt), stat = "identity") +
  coord_flip() +
  labs(y = "type of weighting", 
       x = "Incidence rates per 1,000 person-years") + 
  theme_bw() + 
  theme(legend.position = "none") +
  scale_alpha_discrete(range = c(0.5, 1)) +
  facet_wrap( ~ H_ethn + age_adj)

# save the plot
ggsave(
  filename = paste0(
    path_to_box, path_to_output, "figures/linear vs bspline age/", 
    "IR_x_ethn_linear_vs_bspline_age.png"),
  plot = last_plot()
  # height = 8, units = "in"
)
