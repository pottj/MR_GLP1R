#' ---
#' title: "MR analysis"
#' subtitle: ""
#' author: "Harshika Mohan Raj"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
#' # Load data ####
#' ***

# install.packages("remotes") # Run if remotes package not installed
library(remotes)
install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)
library(tidyverse)    # Data wrangling 
library(TwoSampleMR)  # MR 
library(gt)
library(qqman)

#visual manhattanplot
#manhattan(BMI_fem, chr = "CHR", bp = "POS", p = "P", snp = "SNP", suggestiveline = -log10(1e-04), logp = TRUE)
#BMI_fem %>% filter(P < 1e-04)
#manhattan(PCOS, chr = "chromosome", bp = "base_pair_location", p = "effect_allele_frequency", snp = "effect_allele", suggestiveline = -log10(1e-04), logp = TRUE)

# 1 - DIRECT CALCULATIONS 
#Wald ratio - corresponds to the log odds ratio for the outcome per unit change of the exposure.
#ratio of coefficients, or the Wald ratio -estimating the causal effect of the exposure on the outcome
wald_ratio <- outcome_data$beta/exposure_data$BETA
wald_ratio_standard_error <- outcome_data$standard_error/exposure_data$BETA
z_statistic <- wald_ratio/wald_ratio_standard_error
p_value <- 2*pnorm(abs(z_statistic) ,lower.tail=F)
x<- wald_ratio
y<- wald_ratio_standard_error
z_statistic
p_value

#IVW
IVW_weights <- outcome_data$standard_error^(-2)
inverse_weighted_LR <- lm(outcome_data$beta ~ exposure_data$BETA
                          - 1 ,weights=IVW_weights)
summary(inverse_weighted_LR) # intercept term here is zero in order to calculate the IVW estimate

#MR-Egger Regression
MR_egger_regression <- lm(outcome_data$beta ~ exposure_data$BETA,
                          weights=1/IVW_weights)
summary(MR_egger_regression) 

#IVW alt
bx <- BMI_fem$BETA
bxse <- BMI_fem$SE
by <- PCOS$beta
byse <- PCOS$standard_error
sum(by*bx*byse^-2)/sum(bx^2*byse^-2)


# 2 - MR PACKAGE METHOD 
#BMI female (exposure) vs PCOS (outcome)
filt_fem = BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290")
filt_PCOS = PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290")
mr_obj_fem = mr_input(bx = as.numeric(BMI_fem$BETA[filt_fem]),
                      bxse = as.numeric(BMI_fem$SE[filt_fem]),
                      by = as.numeric(PCOS$beta[filt_PCOS]),
                      byse = as.numeric(PCOS$standard_error[filt_PCOS]),
                      )
       
mr_ivw(mr_input(bx = as.numeric(BMI_fem$BETA[filt_fem]),bxse = as.numeric(BMI_fem$SE[filt_fem]),by = as.numeric(PCOS$beta[filt_PCOS]), byse = as.numeric(PCOS$standard_error[filt_PCOS]))
mr_allmethods(mr_obj_fem)
mr_plot(mr_obj_fem)
mr_loo(mr_obj_fem)
mr_forest(mr_obj_fem, ordered=TRUE)

#BMI male (exposure) vs PCOS (outcome)
filt_mal = BMI_mal$rsID %in% c("rs17757975","rs2254336","rs4714290")
filt_PCOS = PCOS$rs_id %in% c("rs17757975","rs2254336","rs4714290")
mr_obj_mal = mr_input(bx = as.numeric(BMI_mal$BETA[filt_mal]),
                      bxse = as.numeric(BMI_mal$SE[filt_mal]),
                      by = as.numeric(PCOS$beta[filt_PCOS]),
                      byse = as.numeric(PCOS$standard_error[filt_PCOS]),
                      )

#mr_ivw(mr_obj_mal)
mr_allmethods(mr_obj_mal)
mr_plot(mr_obj_mal)
mr_loo(mr_obj_mal)
mr_forest(mr_obj_mal, ordered=TRUE)   

#BMI combined (exposure) vs PCOS (outcome)
filt_com = BMI_comb$rsID %in% c("rs17757975", "rs4714290")
filt_PCOS = PCOS$rs_id %in% c("rs17757975","rs4714290")
mr_obj_comb = mr_input(bx = as.numeric(BMI_comb$BETA[filt_com]),
                      bxse = as.numeric(BMI_comb$SE[filt_com]),
                      by = as.numeric(PCOS$beta[filt_PCOS]),
                      byse = as.numeric(PCOS$standard_error[filt_PCOS]),
                      )

#mr_ivw(mr_obj_comb)
mr_allmethods(mr_obj_comb)
mr_plot(mr_obj_comb)
mr_loo(mr_obj_comb)
mr_forest(mr_obj_comb, ordered=TRUE)   


#' # Save ####
#' ***
#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
