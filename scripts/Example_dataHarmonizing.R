#' ---
#' title: "Example script for data harmonization"
#' subtitle: ""
#' author: "Janne Pott and Harshika Mohan Raj (edits)"
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
#' In this example script, I load data for an exposure (BMI) and outcome (PCOS), and harmonize the effect alleles. 
#' 
#' All necessary R packages and paths to data files are in the source file
#' 
#' This is an classic R file that can be compiled into a report just like an Rmarkdown file. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

# GWAS summary statistics for BMI

Pulit_BMI_female = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit_BMI_male = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit_BMI_combined = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt"

# GWAS summary statistics for PCOS

Venkatesh_PCOS = "/Users/harshikamohanraj/Downloads/GCST90483500.h.tsv"


#' # Load data ####
#' ***
BMI_fem = fread( "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.females.23May2018.txt")
BMI_mal = fread("/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.males.23May2018.txt")
BMI_comb = fread("/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt")

PCOS = fread("/Users/harshikamohanraj/Downloads/PCOS_sumstats.tsv")

#' # Filter data ####
#' ***
#' location of GLPR1: 
#' 
#' - chr6:39,016,574-39,055,519 (GRCh37/hg19 by Ensembl)
#' - chr6:39,048,781-39,091,303 (GRCh38/hg38)
#' 
names(BMI_fem)
BMI_fem = BMI_fem[CHR==6,]
BMI_fem = BMI_fem[POS > 39016574 - 1e6,]
BMI_fem = BMI_fem[POS < 39055519 + 1e6,]

names(BMI_mal)
BMI_mal = BMI_mal[CHR==6,]
BMI_mal = BMI_mal[POS > 39016574 - 1e6,]
BMI_mal = BMI_mal[POS < 39055519 + 1e6,]

names(BMI_comb)
BMI_comb = BMI_comb[CHR==6,]
BMI_comb = BMI_comb[POS > 39016574 - 1e6,]
BMI_comb = BMI_comb[POS < 39055519 + 1e6,]

names(PCOS)
PCOS = PCOS[chromosome==6,]
PCOS = PCOS[base_pair_location > 39016574 - 1e6,]
PCOS = PCOS[base_pair_location < 39055519 + 1e6,]

#' Remove variants with very low frequency
BMI_fem = BMI_fem[Freq_Tested_Allele >=0.01 & Freq_Tested_Allele <=0.99,]
BMI_mal = BMI_mal[Freq_Tested_Allele >=0.01 & Freq_Tested_Allele <=0.99,]
BMI_comb = BMI_comb[Freq_Tested_Allele >=0.01 & Freq_Tested_Allele <=0.99,]
PCOS = PCOS[effect_allele_frequency >=0.01 & effect_allele_frequency <=0.99,]

#' Check overlap of rsID
BMI_fem[,rsID := gsub(":.*","",SNP)]
BMI_mal[,rsID := gsub(":.*","",SNP)]
BMI_comb[,rsID := gsub(":.*","",SNP)]

PCOS = PCOS[rs_id %in% BMI_fem$rsID & rs_id %in% BMI_mal$rsID & rs_id %in% BMI_comb$rsID,]
BMI_fem = BMI_fem[rsID %in% PCOS$rs_id, ]
BMI_mal = BMI_mal[rsID %in% PCOS$rs_id, ]
BMI_comb = BMI_comb[rsID %in% PCOS$rs_id, ]

# Check duplicates & triallelic SNPs
BMI_fem[duplicated(SNP),]
BMI_fem[duplicated(rsID),]
BMI_fem[duplicated(POS),]

#' Check order of data sets
setorder(PCOS,base_pair_location) #Alternative: match command - using id variable to order the other datasets 
stopifnot(BMI_comb$rsID == PCOS$rs_id)
stopifnot(BMI_comb$rsID == BMI_fem$rsid)
stopifnot(BMI_comb$rsID == BMI_mal$rsid)

#' **Summary**: There are 6125 SNPs at the _GLP1R_ locus that are 
#' 
#' - available in females, males, and sex-combined BMI data and PCOS,
#' - have MAF>0.01
#' - are not triallelic or duplicated
#'
#' # Harmonize alleles ####
#' ***
#' At the moment, the SNPs are just the same rsIDs. But I want the same effect allele throughout. 
table(BMI_comb$Tested_Allele == BMI_fem$Tested_Allele,
      BMI_comb$Other_Allele == BMI_fem$Other_Allele)
table(BMI_comb$Tested_Allele == BMI_mal$Tested_Allele,
      BMI_comb$Other_Allele == BMI_mal$Other_Allele)

    #EA and OA -> renaming the columns 

plot(BMI_fem$Freq_Tested_Allele, BMI_mal$Freq_Tested_Allele)
plot(BMI_fem$Freq_Tested_Allele, BMI_comb$Freq_Tested_Allele)

table(BMI_comb$Tested_Allele == PCOS$effect_allele,
      BMI_comb$Other_Allele == PCOS$other_allele)
table(BMI_comb$Tested_Allele == PCOS$other_allele,
      BMI_comb$Other_Allele == PCOS$effect_allele)

plot(BMI_comb$Freq_Tested_Allele, PCOS$effect_allele_frequency)

#' In the PCOS data, there are 3,383 SNPs with same allele coding as the BMI data, and 2,742 SNPs with switch alleles. These SNPs will be flipped back to the BMI coding. 
#' 
filt = BMI_comb$Tested_Allele == PCOS$other_allele & 
  BMI_comb$Other_Allele == PCOS$effect_allele
table(filt)
PCOS[filt,beta := beta * (-1)]
PCOS[filt,effect_allele_frequency := 1-effect_allele_frequency]
PCOS[filt,effect_allele := BMI_comb[filt,Tested_Allele]]
PCOS[filt,other_allele := BMI_comb[filt,Other_Allele]]

#' Now check the transformation with the EAF plot again. 
plot(BMI_comb$Freq_Tested_Allele, PCOS$effect_allele_frequency)

#' Remove the outlier that have different allele frquencies. 
filt = abs(BMI_fem$Freq_Tested_Allele - PCOS$effect_allele_frequency) >0.1
table(filt)

PCOS = PCOS[!filt,]
BMI_comb = BMI_comb[!filt,]
BMI_mal = BMI_mal[!filt,]
BMI_fem = BMI_fem[!filt,]

#' # Save ####
#' ***
#' Save the harmonized data 
#' 
save(BMI_comb,BMI_fem,BMI_mal,PCOS,file = "/Users/harshikamohanraj/Downloads/Input_harmonized.RData")



#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

