# This script checks covariate balance between unweighted RPGEH and CHIS,
# and then generates weights for RPGEH and check covariate balance iteratively. 
# Many of the results are exploratory, and are NOT used for estimation of 
# incidence rates.

# set up packages ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("haven", "tidyverse", "twang", "RColorBrewer", "mitools", "ggpubr", 
       "patchwork"
       # "mice", rlang"
       # "magrittr", "foreign", "ggplot2", "table1", "labelled"
       # "survey", "tableone", "openxlsx", "survival", "mgcv", "miceadds"
       # "openxlsx", "lmtest"
)

options(scipen = 999, digits = 8)

# set up paths and data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/"
path_to_datasets <- "Asian_Americans_dementia_data/aa_selection/"
path_to_output <- "Asian_Americans_dementia/Manuscripts/AA_selection/Code/cleaned_scripts/output/"
chis_harm <- readRDS(paste0(path_to_box, path_to_datasets, 
                            "chis_2005_to_2009_harmonized.RDS"))
rpgeh_harm_imp <- readRDS(paste0(path_to_box, path_to_datasets, 
                                 "rpgeh_imp_harmonized.RDS")) 

# concatenate and relabel RPGEH and CHIS ----
glimpse(chis_harm)
glimpse(rpgeh_harm_imp)
# the only difference is the imp variable (indexing #imputation) in rpgeh

source(paste0(path_to_box, 
              "Asian_Americans_dementia/Manuscripts/AA_selection/Code/",
              "cleaned_scripts/functions_summary.R"))

# main analysis will be excluding Medicare only subjects (for now)
chis_rpgeh <- chis_harm %>% 
  mutate(imp = 0) %>% 
  add_row(rpgeh_harm_imp) %>% 
  t1_relabel()

# vars to check balance for
harmonized_norace <- c(
  "H_age", "H_female", "H_usborn", "H_edu_4", "H_work", "H_retired", 
  "H_marit", "H_hhsize_3", "H_income_pp", "H_health_3", "H_bmi", 
  "H_diab", "H_hyp", "H_smk"
)

# calculate unconditional RPGEH sampling odds ----
temp <- chis_rpgeh %>% filter(imp %in% c(0,1))
# overall 
p_rpgeh <- sum(temp$in_rpgeh) / sum(temp$smplwt)
uncon_rpgeh_odds <- p_rpgeh / (1 - p_rpgeh)
# by ethnicity
p_rpgeh_ethn <- temp %>% group_by(H_ethn) %>%
  summarize(marginal_p = sum(in_rpgeh) / sum(smplwt),
            odds = marginal_p / (1 - marginal_p)) 
remove(temp)

# set up function for checking covariate balance and generating weights ----
# given data and a specified ethnicity, 
# (1) if unweighted = TRUE, we check the unweighted covariate balance and 
# save the bal.stat output for the 40 imputations and the std mean diff plot
# (2) if unweighted = FALSE, we also supply a formula to generate weights with.
# In addition to the output in the unweighted case, numerical and graphical 
# summaries of weights and propensity scores are generated using CHIS and the 
# last imputation of RPGEH. 


wt_cov_bal <- function(data, ethn, fml = NULL, unweighted = FALSE) {
  # testing
  # data <- chis_rpgeh
  # ethn <- "Chinese"
  # fml <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
  # unweighted <- FALSE
  
  # subset to a specific ethnicity
  ethn_subset <- data %>% filter(H_ethn == ethn)
  
  if (unweighted) {
    # check covariate balance for unweighted RPGEH
    # set up a matrix to store std.eff.sz from each imputation
    for (i in 1:40) {
      temp <- ethn_subset %>% 
        # select the one copy of CHIS and the i-th imputed RPGEH
        filter(imp %in% c(0, i)) %>% 
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
      
      # store std.eff.sz for each imputation
      if (i == 1) {
        std_eff_sz <- tibble(
          var = rownames(temp$results),
          imp_1 = -temp$results$std.eff.sz # reverse the sign
        )
      } else {
        # reverse the sign
        std_eff_sz[, paste0("imp_", i)] <- -temp$results$std.eff.sz
      }
    }
  } else {
    # if needs weighting, fit glm and generate weights, etc
    
    ## fit glm model for each imputation given formula
    for (i in 1:40) {
      mod <- ethn_subset %>% filter(imp %in% c(0, i)) %>%
        glm(
          as.formula(fml),
          data = .,
          family = binomial(link = "logit"),
          weights = .$smplwt
        )
      ## calculate propensity scores
      ethn_subset[which(ethn_subset$imp %in% c(0, i)), "p1"] <- 
        predict.glm(mod, type = "response")
    }
    # note at at this point, ethn_subset has ps for each imputation of RPGEH
    # while for CHIS only the ps from the glm model with the last imputation 
    # of RPGEH are stored
    
    ## calculate weights
    ## pull the ethnicity specific RPGEH odds
    rpgeh_odds <- p_rpgeh_ethn %>% filter(H_ethn == ethn) %>% pull(odds)
    ethn_subset <- ethn_subset %>% mutate(
      iow1 = (1 - p1) / p1,
      # note that for CHIS, sw1 = 1 and sw1_twang = chis survey weight
      sw1 = ifelse(in_rpgeh == 1, iow1 * rpgeh_odds, 1),
      sw1_twang = sw1 * smplwt
    )
    
    ## summary stats for stabilized weights looking at last copy of RPGEH
    ## used to label the density plots
    sw.sum <- summary(ethn_subset[which(ethn_subset$imp == 40),]$sw1)
    sw.sum["sd"] <- sd(ethn_subset[which(ethn_subset$imp == 40),]$sw1)
    sw.sum <- round(sw.sum, 3)
    plt_labels <- c(
      paste0(sw.sum[4], " (", sw.sum[7], ")"), # mean (sd)
      # 1st and 3rd quartiles
      paste0(sw.sum[2]), paste0(sw.sum[5])  
    )
    
    # density plots for propensity scores and weights 
    ## for CHIS and the last copy of RPGEH
    plot_data <- ethn_subset %>% filter(imp %in% c(0, 40))
    ## propensity scores for both studies
    plt_ps <- ggplot(plot_data, aes(x = p1, group = in_chis)) +
      geom_density(aes(color = factor(in_chis))) +
      labs(x = "Propensity Score", color = "In CHIS")
    ## weights for RPGEH
    plt_wt <- plot_data %>% filter(in_chis == 0) %>% 
      ggplot(., aes(x = sw1, group = in_chis)) +
      geom_density(aes(color = factor(in_chis))) +
      # add summary statistics as vertical lines and labels
      geom_vline(
        xintercept = c(sw.sum[4], sw.sum[2], sw.sum[5]),
        linetype = c("solid", "dashed", "dashed"), color = "grey") + 
      annotate(geom = 'text', 
               label = plt_labels, 
               x = c(sw.sum[4], sw.sum[2], sw.sum[5]), 
               y = c(Inf, 0, 0), 
               hjust = 0, vjust = c(2, 0, 0), size = 3) +
      labs(x = "Stabilized Weight", color = "In CHIS")
    ## combine the two plots
    ps.wt.plot <- ggarrange(plt_ps, plt_wt, nrow = 2, 
                            common.legend = TRUE, legend = "bottom")
    
    # check covariate balance for weighted RPGEH
    # set up a matrix to store std.eff.sz from each imputation
    for (i in 1:40) {
      temp <- ethn_subset %>% 
        filter(imp %in% c(0, i)) %>% 
        as.data.frame() %>% 
        bal.stat(
          data = .,
          vars = harmonized_norace,
          w.all = .[, "sw1_twang"],
          treat.var = 'in_chis',
          sampw = .$smplwt,
          estimand = "ATT",
          multinom = F
        )
      
      # store std.eff.sz for each imputation
      if (i == 1) {
        std_eff_sz <- tibble(
          var = rownames(temp$results),
          imp_1 = -temp$results$std.eff.sz # reverse the sign
        )
      } else {
        # reverse the sign
        std_eff_sz[, paste0("imp_", i)] <- -temp$results$std.eff.sz
      }
    }
  }
  
  # calculate mean std effect size across imputations
  std_eff_sz[, "mean"] <- rowMeans(std_eff_sz[, 2:41])
  # plot std eff sz 
  cov.bal.plot <- std_eff_sz %>% 
    filter(
      # remove some redundant rows for simplicity
      str_detect(var, "No|Male", negate = TRUE)
    ) %>% 
    ggplot(., aes(x = var, y = mean)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    labs(y = "Standardized mean difference\n(RPGEH-CHIS)/SD[CHIS]", 
         x = NULL) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 10),
      legend.position = "right", 
      aspect.ratio = 5 / 3
    ) +
    ylim(-1, 1.25) +
    coord_flip()
  
  if (unweighted) {
    return(list(std_eff_sz = std_eff_sz, cov_bal_plot = cov.bal.plot))
  } else {
    return(list(wt_df = ethn_subset, 
                sw_sum = sw.sum, 
                std_eff_sz = std_eff_sz, 
                # plot_data = plot_data, 
                cov_bal_plot = cov.bal.plot, 
                ps_wt_plot = ps.wt.plot))
  }
}


# Chinese ----
# unweighted
ch0_res <- wt_cov_bal(chis_rpgeh, "Chinese", unweighted = TRUE)

# weighted with basic demographic variables
# common first step for all ethnicities 
fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
ch1_res <- wt_cov_bal(chis_rpgeh, "Chinese", fml1)
# HHsize and health are not balanced: added to fml2

fml2 <- paste0(fml1, " + H_hhsize_3 + H_health_3")
ch2_res <- wt_cov_bal(chis_rpgeh, "Chinese", fml2)
# hypertension is not balanced: added to fml3

fml3 <- paste0(fml2, " + H_hyp")
ch3_res <- wt_cov_bal(chis_rpgeh, "Chinese", fml3)
# this is the final model

# nonlinear terms in age are explored, but no significant improvement
# code omitted

# save the output objects
# We only save one unweighted and one weighted function output here; 
# the weighted one uses the final glm formula that we decided on. 
# This is used later to generate figures for the manuscript. 

# save(ch0_res, ch3_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/ch_res.RData"))

# save the plots for comparisons:
# these plots are only for exploration and discussion, not used for manuscript; 
# they compare the covariate balance plots and density plots for ps and weights 
# side by side across different models

# all models
# mods <- paste0("Chinese", "\n",
#                "Level 1: ",
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ",
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 3: ",
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# ggsave(
#   plot = ggarrange(ch0_res$cov_bal_plot,
#                    ch1_res$cov_bal_plot,
#                    ch2_res$cov_bal_plot,
#                    ch3_res$cov_bal_plot,
#                    ggplot() + theme_void() +
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    ch1_res$ps_wt_plot,
#                    ch2_res$ps_wt_plot,
#                    ch3_res$ps_wt_plot,
#                    nrow = 2, ncol = 4),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_ch.png"),
#   width = 20, height = 8, units = "in", bg = "white"
# )

# Filipino ----
fi0_res <- wt_cov_bal(chis_rpgeh, "Filipino", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
fi1_res <- wt_cov_bal(chis_rpgeh, "Filipino", fml1)

# add HHsize since it's not balanced
fml2 <- paste0(fml1, " + H_hhsize_3")
fi2_res <- wt_cov_bal(chis_rpgeh, "Filipino", fml2)

# add hypertension and bmi since these are still not balanced
fml3 <- paste0(fml2, " + H_hyp + H_bmi")
fi3_res <- wt_cov_bal(chis_rpgeh, "Filipino", fml3)

# add retired to further improve balance
fml4 <- paste0(fml3, " + H_retired")
fi4_res <- wt_cov_bal(chis_rpgeh, "Filipino", fml4)

# add health since it is an important covariate for dementia
fml5 <- paste0(fml4, " + H_health_3")
fi5_res <- wt_cov_bal(chis_rpgeh, "Filipino", fml5)
# this is the final model

# save the output objects
# save(fi0_res, fi5_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/fi_res.RData"))

# save the plots
# mods <- paste0("Filipino", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n", 
#                "Level 3: ", 
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n", 
#                "Level 4: ", 
#                substring(fml4, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# 
# ggsave(
#   plot = ggarrange(fi0_res$cov_bal_plot, 
#                    fi1_res$cov_bal_plot, 
#                    fi2_res$cov_bal_plot, 
#                    fi3_res$cov_bal_plot, 
#                    fi4_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    fi1_res$ps_wt_plot, 
#                    fi2_res$ps_wt_plot,
#                    fi3_res$ps_wt_plot,
#                    fi4_res$ps_wt_plot,
#                    nrow = 2, ncol = 5),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_fi.png"),
#   width = 25, height = 8, units = "in", bg = "white"
# )

# Japanese ----
ja0_res <- wt_cov_bal(chis_rpgeh, "Japanese", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
ja1_res <- wt_cov_bal(chis_rpgeh, "Japanese", fml1)

# add HHsize since it's not balanced
fml2 <- paste0(fml1, " + H_hhsize_3")
ja2_res <- wt_cov_bal(chis_rpgeh, "Japanese", fml2)

# add hypertension since it's not balanced
# add health since it's generally important
fml3 <- paste0(fml2, " + H_hyp + H_health_3")
ja3_res <- wt_cov_bal(chis_rpgeh, "Japanese", fml3)
# this is the final model

# save the output objects
# save(ja0_res, ja3_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/ja_res.RData"))

# save the plots
# mods <- paste0("Japanese", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_"), "\n", 
#                "Level 3: ", 
#                substring(fml3, 12) %>% str_remove_all("H_"))
# 
# ggsave(
#   plot = ggarrange(ja0_res$cov_bal_plot, 
#                    ja1_res$cov_bal_plot, 
#                    ja2_res$cov_bal_plot, 
#                    ja3_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    ja1_res$ps_wt_plot, 
#                    ja2_res$ps_wt_plot, 
#                    ja3_res$ps_wt_plot, 
#                    nrow = 2, ncol = 4),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_ja.png"),
#   width = 20, height = 8, units = "in", bg = "white"
# )

# Korean ----
ko0_res <- wt_cov_bal(chis_rpgeh, "Korean", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
ko1_res <- wt_cov_bal(chis_rpgeh, "Korean", fml1)

# HHsize and health are not balanced; adding here
fml2 <- paste0(fml1, " + H_hhsize_3 + H_health_3")
ko2_res <- wt_cov_bal(chis_rpgeh, "Korean", fml2)

# adding both HHsize and health appears to worsen the balance of other vars
# try adding these two individually:
fml2.1 <- paste0(fml1, " + H_hhsize_3")
ko2.1_res <- wt_cov_bal(chis_rpgeh, "Korean", fml2.1)

fml2.2 <- paste0(fml1, " + H_health_3")
ko2.2_res <- wt_cov_bal(chis_rpgeh, "Korean", fml2.2)

# adding hhsize and health individually does not affect other covariates
# adding them both at the same time does

# continue on by adding hypertension, which is unbalanced
# fml3: add hyp to model with both HHsize and health
fml3 <- paste0(fml2, " + H_hyp")
ko3_res <- wt_cov_bal(chis_rpgeh, "Korean", fml3)

# fml3.1 and 3.2: add hyp to model with either HHsize or health
fml3.1 <- paste0(fml2.1, " + H_hyp")
ko3.1_res <- wt_cov_bal(chis_rpgeh, "Korean", fml3.1)

fml3.2 <- paste0(fml2.2, " + H_hyp")
ko3.2_res <- wt_cov_bal(chis_rpgeh, "Korean", fml3.2)

# since we believe health is more important than HHsize for dementia, 
# the final model decided on is fml3.2

# save the output objects
# save(ko0_res, ko3.2_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/ko_res.RData"))

# save the plots
# mods <- paste0("Korean", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 3: ", 
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# 
# ggsave(
#   plot = ggarrange(ko0_res$cov_bal_plot, 
#                    ko1_res$cov_bal_plot, 
#                    ko2_res$cov_bal_plot, 
#                    ko3_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    ko1_res$ps_wt_plot, 
#                    ko2_res$ps_wt_plot, 
#                    ko3_res$ps_wt_plot, 
#                    nrow = 2, ncol = 4),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_ko.png"),
#   width = 20, height = 8, units = "in", bg = "white"
# )
# 
# text1 <- paste0("Column 1: demographics\n", 
#                 "Column 2: add hhsize (top) and health (bottom) individually")
# text2 <- paste0("Column 3: add hhsize and health together\n", 
#                 "Column 4: add hyp to hhsize (top) and health (bottom) individually\n", 
#                 "Column 5: add hyp to hhsize and health together")
# ggsave(
#   plot = ggarrange(ko1_res$cov_bal_plot, 
#           ko2.1_res$cov_bal_plot, 
#           ko2_res$cov_bal_plot, 
#           ko3.1_res$cov_bal_plot,
#           ko3_res$cov_bal_plot,
#           ggplot() + theme_void() + 
#             annotate("text", x = 1, y = 1, label = text1, size = 4), 
#           ko2.2_res$cov_bal_plot,
#           ggplot() + theme_void() + 
#             annotate("text", x = 1, y = 1, label = text2, size = 4),
#           ko3.2_res$cov_bal_plot, 
#           nrow = 2, ncol = 5), 
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_ko_add.png"),
#   width = 25, height = 8, units = "in", bg = "white"
# )


# Pacific Islander ----
pi0_res <- wt_cov_bal(chis_rpgeh, "Pacific Islander", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
pi1_res <- wt_cov_bal(chis_rpgeh, "Pacific Islander", fml1)

# health and hypertension are still not balanced; added
fml2 <- paste0(fml1, " + H_health_3 + H_hyp")
pi2_res <- wt_cov_bal(chis_rpgeh, "Pacific Islander", fml2)

# although smk is not balanced, adding it does not improve its balance or 
# the balance of other covariates
fml3 <- paste0(fml2, " + H_smk")
pi3_res <- wt_cov_bal(chis_rpgeh, "Pacific Islander", fml3)

# the final model is fml2

# save the output objects
# save(pi0_res, pi2_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/pi_res.RData"))

# save the plots
# mods <- paste0("Pacific Islander", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 3: ", 
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# 
# ggsave(
#   plot = ggarrange(pi0_res$cov_bal_plot, 
#                    pi1_res$cov_bal_plot, 
#                    pi2_res$cov_bal_plot, 
#                    pi3_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    pi1_res$ps_wt_plot, 
#                    pi2_res$ps_wt_plot, 
#                    pi3_res$ps_wt_plot, 
#                    nrow = 2, ncol = 4),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_pi.png"),
#   width = 20, height = 8, units = "in", bg = "white"
# )

# South Asian ----
sa0_res <- wt_cov_bal(chis_rpgeh, "South Asian", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
sa1_res <- wt_cov_bal(chis_rpgeh, "South Asian", fml1)

# retired and health are not balanced; added here
fml2 <- paste0(fml1, " + H_retired + H_health_3")
sa2_res <- wt_cov_bal(chis_rpgeh, "South Asian", fml2)

# hypertension is not balanced; added here
fml3 <- paste0(fml2, " + H_hyp")
sa3_res <- wt_cov_bal(chis_rpgeh, "South Asian", fml3)

# adding HHsize does not improve overall balance
fml4 <- paste0(fml3, " + H_hhsize_3")
sa4_res <- wt_cov_bal(chis_rpgeh, "South Asian", fml4)

# when running bootstraps, the small number of usborn South Asians 
# are causing computation problems -- decided to drop it from the final model
fml5 <- "in_rpgeh ~ H_age + H_female + H_edu_4 + H_retired + H_health_3 + H_hyp"
sa5_res <- wt_cov_bal(chis_rpgeh, "South Asian", fml5)

# save the output objects
# save(sa0_res, sa5_res,
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/sa_res.RData"))

# save the plots
# mods <- paste0("South Asian", "\n",
#                "Level 1: ",
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ",
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 3: ",
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 4: ",
#                substring(fml4, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 5: ",
#                substring(fml5, 12) %>% str_remove_all("H_") %>% str_wrap(29)
#                )
# 
# ggsave(
#   plot = ggarrange(sa0_res$cov_bal_plot,
#                    sa1_res$cov_bal_plot,
#                    sa2_res$cov_bal_plot,
#                    sa3_res$cov_bal_plot,
#                    sa4_res$cov_bal_plot,
#                    sa5_res$cov_bal_plot, 
#                    ggplot() + theme_void() +
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    sa1_res$ps_wt_plot,
#                    sa2_res$ps_wt_plot,
#                    sa3_res$ps_wt_plot,
#                    sa4_res$ps_wt_plot,
#                    sa5_res$ps_wt_plot,
#                    nrow = 2, ncol = 6),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_sa.png"),
#   width = 30, height = 8, units = "in", bg = "white"
# )

# Vietnamese ----
vi0_res <- wt_cov_bal(chis_rpgeh, "Vietnamese", unweighted = TRUE)
# usborn is missing since no one was US born
# fair or poor health is out of bounds, at ~ -1.1

# fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
fml1 <- "in_rpgeh ~ H_age + H_female + H_edu_4" # drop usborn (09/30)
# usborn is dropped (09/30) since no one was US born in either study
vi1_res <- wt_cov_bal(chis_rpgeh, "Vietnamese", fml1)

# work and health are unbalanced; added
fml2 <- paste0(fml1, " + H_work + H_health_3")
vi2_res <- wt_cov_bal(chis_rpgeh, "Vietnamese", fml2)

# adding smk to level 3 only makes smk more unbalanced; not done
# bmi is the next most unbalanced variable; added
fml3 <- paste0(fml2, " + H_bmi")
vi3_res <- wt_cov_bal(chis_rpgeh, "Vietnamese", fml3)

# adding hypertension last
fml4 <- paste0(fml3, " + H_hyp")
vi4_res <- wt_cov_bal(chis_rpgeh, "Vietnamese", fml4)
# final model

# save the output objects
# save(vi0_res, vi4_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/vi_res.RData"))

# load(paste0(path_to_box, path_to_datasets, "weighted_data/vi_res.RData"))

# save the plots
# mods <- paste0("Vietnamese", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n",
#                "Level 3: ", 
#                substring(fml3, 12) %>% str_remove_all("H_") %>% str_wrap(29), "\n", 
#                "Level 4: ", 
#                substring(fml4, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# 
# ggsave(
#   plot = ggarrange(vi0_res$cov_bal_plot, 
#                    vi1_res$cov_bal_plot, 
#                    vi2_res$cov_bal_plot, 
#                    vi3_res$cov_bal_plot, 
#                    vi4_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    vi1_res$ps_wt_plot, 
#                    vi2_res$ps_wt_plot, 
#                    vi3_res$ps_wt_plot, 
#                    vi4_res$ps_wt_plot, 
#                    nrow = 2, ncol = 5),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_vi.png"),
#   width = 25, height = 8, units = "in", bg = "white"
# )

# Other SE Asian ----
# small sample size - dropped from analysis

ose0_res <- wt_cov_bal(chis_rpgeh, "Other SE Asian", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
ose1_res <- wt_cov_bal(chis_rpgeh, "Other SE Asian", fml1)

# want to add H_hhsize_3, but glm gives warnings
fml2 <- paste0(fml1, " + H_hyp")
ose2_res <- wt_cov_bal(chis_rpgeh, "Other SE Asian", fml2)

# save the output objects
# save(ose0_res, ose2_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/ose_res.RData"))

# save the plots
# mods <- paste0("Other SE Asian", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n",
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_") %>% str_wrap(29))
# 
# ggsave(
#   plot = ggarrange(ose0_res$cov_bal_plot, 
#                    ose1_res$cov_bal_plot, 
#                    ose2_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    ose1_res$ps_wt_plot, 
#                    ose2_res$ps_wt_plot,
#                    nrow = 2, ncol = 3),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_ose.png"),
#   width = 15, height = 8, units = "in", bg = "white"
# )

# Multiple Ethnicities ----
# small sample size - dropped from analysis

mu0_res <- wt_cov_bal(chis_rpgeh, "Multiple Ethnicities", unweighted = TRUE)

fml1 <- "in_rpgeh ~ H_age + H_female + H_usborn + H_edu_4"
mu1_res <- wt_cov_bal(chis_rpgeh, "Multiple Ethnicities", fml1)

fml2 <- paste0(fml1, " + H_smk + H_health_3")
mu2_res <- wt_cov_bal(chis_rpgeh, "Multiple Ethnicities", fml2)

# save the output objects
# save(mu0_res, mu2_res, 
#      file = paste0(path_to_box, path_to_datasets, "weighted_data/mu_res.RData"))

# save the plots
# mods <- paste0("Multiple ethnicities", "\n",
#                "Level 1: ", 
#                substring(fml1, 12) %>% str_remove_all("H_"), "\n", 
#                "Level 2: ", 
#                substring(fml2, 12) %>% str_remove_all("H_"))
# 
# ggsave(
#   plot = ggarrange(mu0_res$cov_bal_plot, 
#                    mu1_res$cov_bal_plot, 
#                    mu2_res$cov_bal_plot, 
#                    ggplot() + theme_void() + 
#                      annotate("text", x = 1, y = 1, label = mods, size = 4),
#                    mu1_res$ps_wt_plot, 
#                    mu2_res$ps_wt_plot, 
#                    nrow = 2, ncol = 3),
#   file = paste0(path_to_box, path_to_output, "figures/cov_bal_mu.png"),
#   width = 15, height = 8, units = "in", bg = "white"
# )





