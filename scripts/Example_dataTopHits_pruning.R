#' ---
#' title: "Example script for top SNP identification"
#' subtitle: ""
#' author: "Janne Pott and Harshika Mohan Raj"
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
#' In this example script, I load data the harmonized data for BMI and check for potential instruments in the three data sets. 
#' 
#' All necessary R packages and paths to data files are in the source file
#' 
#' This is an classic R file that can be compiled into a report just like an Rmarkdown file. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
#' # Load data ####
#' ***

#' # Significant SNPs ####
#' ***
#' Check number of significant SNPs
BMI_fem[,table(P<5e-8,P<1e-6)]
BMI_mal[,table(P<5e-8,P<1e-6)]
BMI_comb[,table(P<5e-8,P<1e-6)]

#' Add setting parameter for later re-identification
BMI_fem[,setting := "females"]
BMI_mal[,setting := "males"]
BMI_comb[,setting := "combined"]

#' Combine all settings, and filter for genome-wide significant SNPs
BMI = rbind(BMI_fem,BMI_mal, BMI_comb)
BMI = BMI[P<5e-8,]
BMI[,length(unique(rsID))]
BMI[,table(setting)]

#' There are 42 unique SNPs that reach genome-wide significance in at least one setting. 
#' 
#' # Pruning ####
#' ***
#' Now I print the rsIDs of the 42 SNPs and copy them into LDlink: https://ldlink.nih.gov/?tab=ldmatrix - using European as reference population. 
#' 
BMI$rsID

#' There are four clusters, nicely sorted by position. Hence, I add a cluster variable by position information
setorder(BMI,POS)
BMI[,cluster := "cluster1"]
BMI[POS<=38227592,cluster := "cluster1"]
BMI[POS>38227592 & POS<= 39044558,cluster := "cluster2"]
BMI[POS>39044558 & POS<= 39920827,cluster := "cluster3"]
BMI[POS>39920827,cluster := "cluster4"]

#' Now I get the best SNP per cluster (lowest p-value)
topList = copy(BMI)
setorder(topList,P)
topList = topList[!duplicated(cluster)]
topList[,table(setting)]
setorder(topList,POS)
topList

#' # Check per setting ####
#' ***
BMI_comb[rsID %in% topList$rsID]
BMI_mal[rsID %in% topList$rsID]
BMI_fem[rsID %in% topList$rsID]

#' **Summary**:
#' 
#' - the two hits from the combined data are not genome-wide significant in the sex-stratified data sets, probably due to lower sample size = lower power; the effect sizes are similar
#' - the hit from the sex-stratified data is sex-specific, e.g. the effect sizes differ between the sexes and in the combined setting the effect is averaged over both (higher variability and greater SE). 
#' 
#' # Save ####
#' ***
#' Save the harmonized data 
#' 
save(BMI_comb,BMI_fem,BMI_mal,PCOS,file = "/Users/harshikamohanraj/Downloads/Input_TopHitsPruned.RData")

names(BMI_comb)
# install.packages("remotes") # Run if remotes package not installed
library(remotes)
install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)

library(tidyverse)    # Data wrangling 
library(TwoSampleMR)  # MR 
library(gt)

#' MR analysis 
PCOS
exposure_data <- BMI_fem
outcome_data <- PCOS

library(qqman)
manhattan(BMI_fem, chr = "CHR", bp = "POS", p = "P", snp = "SNP", suggestiveline = -log10(1e-04), logp = TRUE)
BMI_fem %>% filter(P < 1e-04)

manhattan(PCOS, chr = "chromosome", bp = "base_pair_location", p = "effect_allele_frequency", snp = "effect_allele", suggestiveline = -log10(1e-04), logp = TRUE)

#Wald ratio corresponds to the log odds ratio for the outcome per unit change of the exposure.
#ratio of coefficients, or the Wald ratio -estimating the causal effect of the exposure on the outcome
wald_ratio <- outcome_data$beta/exposure_data$BETA
wald_ratio_standard_error <- outcome_data$standard_error/exposure_data$BETA
z_statistic <- wald_ratio/wald_ratio_standard_error
p_value <- 2*pnorm(abs(z_statistic) ,lower.tail=F)

x<- wald_ratio
y<- wald_ratio_standard_error
z_statistic
p_value

IVW_weights <- outcome_data$standard_error^(-2)
inverse_weighted_LR <- lm(outcome_data$beta ~ exposure_data$BETA
                          - 1 ,weights=IVW_weights)
summary(inverse_weighted_LR) # intercept term here is zero in order to calculate the IVW estimate

#MR-Egger Regression
MR_egger_regression <- lm(outcome_data$beta ~ exposure_data$BETA,
                          weights=1/IVW_weights)
summary(MR_egger_regression) 

#alternative bmi fem 
mr_obj_bf = mr_input(bx = BMI_fem$BETA,
                  bxse = BMI_fem$SE,
                  by = PCOS$beta,
                  byse = PCOS$standard_error,
                  )
#IVW estimate can also be obtained by direct calculation
bx <- BMI_fem$BETA
bxse <- BMI_fem$SE
by <- PCOS$beta
byse <- PCOS$standard_error
sum(by*bx*byse^-2)/sum(bx^2*byse^-2)


mr_allmethods(mr_obj_bf)
# scatter plot
mr_plot(mr_obj_bf)
# single snp plot
mr_forest(mr_obj_bf, ordered=TRUE)
# leave one out plot
mr_loo(mr_obj_bf)
# funnel plot
mr_funnel(mr_obj_bf)

#alternative bmi male 
bmi_snps <- c("rs17757975", "rs116208210", "rs4714290")
mr_obj_bm = mr_input(bx = as.numeric(BMI_mal$BETA),
                     bxse = as.numeric(BMI_mal$SE),
                     by = as.numeric(PCOS$beta),
                     byse = as.numeric(PCOS$standard_error)
                     )
mr_ivw(mr_input(as.numeric(BMI_mal$BETA), as.numeric(BMI_mal$SE), as.numeric(PCOS$beta), as.numeric(PCOS$standard_error))
mr_method_list()
mr_allmethods(mr_obj_bm)
mr_plot(mr_obj_bm, interactive=FALSE, orientate=TRUE)
# single snp plot
mr_forest(mr_obj_bm, ordered=TRUE)
# leave one out plot
mr_loo(mr_obj_bm)
# funnel plot
mr_funnel(mr_obj_bm)

#alternative combined 
BMI_fem[BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290"),]
PCOS[PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290"),]
mr_obj_fem = mr_input(bx = as.numeric(BMI_fem[BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290"),]$BETA),
                     bxse = as.numeric(BMI_fem[BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290"),]$SE),
                     by = as.numeric(PCOS[PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290"),]$beta),
                     byse = as.numeric(PCOS[PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290"),]$standard_error),
)

mr_ivw(mr_input(bx = as.numeric(BMI_fem[BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290"),]$BETA),
                bxse = as.numeric(BMI_fem[BMI_fem$rsID %in% c("rs17757975","rs116208210","rs4714290"),]$SE),
                by = as.numeric(PCOS[PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290"),]$beta),
                byse = as.numeric(PCOS[PCOS$rs_id %in% c("rs17757975","rs116208210","rs4714290"),]$standard_error)))
mr_allmethods(mr_obj_fem)
mr_plot(mr_obj_fem)
mr_loo(mr_obj_fem)
mr_forest(mr_obj_fem, ordered=TRUE)

# scatter plot
snps_of_interest <- c("rs17757975", "rs116208210", "rs4714290")

mr_plot(mr_obj_com, SNP= c("rs17757975", "rs116208210", "rs4714290"))

# single snp plot
mr_forest(mr_obj_com, ordered=TRUE)
# leave one out plot
mr_loo(mr_obj_com)
# funnel plot
mr_funnel(mr_obj_com)

##MR alternative 
harmonised_data = fread("/Users/harshikamohanraj/Downloads/Input_harmonized.RData")
res <- mr(harmonised_data)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


