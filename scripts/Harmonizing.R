#' ---
#' title: "Harmonization"
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
#' Load data for an exposure (BMI) and outcome (PCOS), and harmonize the effect alleles. 
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

Pulit_BMI_female = "/Users/harshikamohanraj/Internship/bmi.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit_BMI_male = "/Users/harshikamohanraj/Internship/bmi.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit_BMI_combined = "/Users/harshikamohanraj/Internship/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt"

# GWAS summary statistics for PCOS
Venkatesh_PCOS = "/Users/harshikamohanraj/Internship/GCST90483500.h.tsv"

#' # Load data ####
#' ***
BMI_fem = fread("/Users/harshikamohanraj/Desktop/Internship/bmi.giant-ukbb.meta-analysis.females.23May2018.txt")
BMI_mal = fread("/Users/harshikamohanraj/Desktop/Internship/bmi.giant-ukbb.meta-analysis.males.23May2018.txt")
BMI_comb = fread("/Users/harshikamohanraj/Desktop/Internship/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt")
PCOS = fread("/Users/harshikamohanraj/Desktop/Internship/PCOS_sumstats.tsv")

#' # Filter data ####
#' ***
#' location of GLPR1: 
#' 
#' - chr6:39,016,574-39,055,519 (GRCh37/hg19 by Ensembl)
#' - chr6:39,048,781-39,091,303 (GRCh38/hg38)
#' 
names(BMI_fem)
head(BMI_fem)
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
save(BMI_comb,BMI_fem,BMI_mal,PCOS,file = "/Users/harshikamohanraj/Desktop/Internship/Input_harmonized.RData")


##Positive controls harmonising - ALL-CAUSE MORTALITY 
all_cause = fread("/Users/harshikamohanraj/Desktop/Internship/lifegen_phase2_bothpl_alldr_2017_09_18.tsv")

all_cause = all_cause[chr==6,]
all_cause = all_cause[pos > 39016574 - 1e6,]
all_cause = all_cause[pos < 39055519 + 1e6,]
#' Remove variants with very low frequency
all_cause = all_cause[freq1 >=0.01 & freq1 <=0.99,]

#' Check overlap of rsID
all_cause = all_cause[rsid %in% BMI_fem$rsID & rsid %in% BMI_mal$rsID & rsid %in% BMI_comb$rsID,]
BMI_fem = BMI_fem[rsID %in% all_cause$rsid, ]
BMI_mal = BMI_mal[rsID %in% all_cause$rsid, ]
BMI_comb = BMI_comb[rsID %in% all_cause$rsid, ]

#' Check order of data sets
setorder(all_cause,pos) #Alternative: match command - using id variable to order the other datasets 
stopifnot(BMI_comb$rsID == all_cause$rsid)
head(all_cause)
#' # Harmonize alleles ####
#' ***
filt = BMI_comb$Tested_Allele == all_cause$a0 & 
  BMI_comb$Other_Allele == all_cause$a1
table(filt)
all_cause[filt,beta1 := beta1 * (-1)]
all_cause[filt,freq1 := 1-freq1]
all_cause[filt,a1 := BMI_comb[filt,Tested_Allele]]
all_cause[filt,a0 := BMI_comb[filt,Other_Allele]]

#' Now check the transformation with the EAF plot again. 
plot(BMI_comb$Freq_Tested_Allele, all_cause$freq1)

head(all_cause)

#' Remove the outlier that have different allele frquencies. 
filt = abs(BMI_fem$Freq_Tested_Allele - all_cause$freq1) >0.1
table(filt)
all_cause = all_cause[!filt,]


##Positive controls harmonising - CAD 
cad = fread("/Users/harshikamohanraj/Desktop/Internship/CAD_GWAS_SEX_STRATIFIED.txt.gz")

head(cad)
head(all_cause)
cad = cad[CHR==6,]
cad = cad[BP > 39016574 - 1e6,]
cad = cad[BP < 39055519 + 1e6,]
#' Remove variants with very low frequency
cad = cad[eaf >=0.01 & eaf <=0.99,]

#' Check overlap of rsID
cad = cad[rsid_ukb %in% BMI_fem$rsID & rsid_ukb %in% BMI_mal$rsID & rsid_ukb %in% BMI_comb$rsID,]
BMI_fem = BMI_fem[rsID %in% cad$rsid_ukb, ]
BMI_mal = BMI_mal[rsID %in% cad$rsid_ukb, ]
BMI_comb = BMI_comb[rsID %in% cad$rsid_ukb, ]

#' Check order of data sets
setorder(cad,BP) #Alternative: match command - using id variable to order the other datasets 
stopifnot(BMI_comb$rsID == cad$rsid_ukb)

#' # Harmonize alleles ####
#' ***
filt = BMI_comb$Tested_Allele == cad$other_allele & 
  BMI_comb$Other_Allele == cad$reference_allele
table(filt)
cad[filt,beta := beta * (-1)]
cad[filt,eaf := 1-eaf]
cad[filt,reference_allele := BMI_comb[filt,Tested_Allele]]
cad[filt,other_allele := BMI_comb[filt,Other_Allele]]

#' Now check the transformation with the EAF plot again. 
plot(BMI_comb$Freq_Tested_Allele, cad$eaf)


#' Remove the outlier that have different allele frquencies. 
filt = abs(BMI_fem$Freq_Tested_Allele - cad$eaf) >0.1
table(filt)
cad = cad[!filt,]
BMI_comb = BMI_comb[!filt,]
BMI_fem = BMI_fem[!filt,]
BMI_mal = BMI_mal[!filt,]

## OUTCOME 2 = HbA1c 
hb = fread("/Users/harshikamohanraj/Desktop/Internship/30750_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz")
PCOS = fread("/Users/harshikamohanraj/Desktop/Internship/PCOS_sumstats.tsv")

#' # Filter data ####
#' ***
#' location of GLPR1: 
#' 
var = unlist(strsplit(hb$variant,":"))
chr = var[seq(1,length(var),4)]
pos = var[seq(2,length(var),4)]
a1 = var[seq(3,length(var),4)]
a2 = var[seq(4,length(var),4)]

hb[,chr := as.numeric(chr)]
hb[,pos := as.numeric(pos)]
hb[,a1 := as.character(a1)]
hb[,a2 := as.character(a2)]
head(hb)

hb = hb[chr== "6",]
hb = hb[pos > 39016574 - 1e6,]
hb = hb[pos < 39055519 + 1e6,]

names(hb)
names(PCOS)
PCOS = PCOS[chromosome==6,]
PCOS = PCOS[base_pair_location > 39016574 - 1e6,]
PCOS = PCOS[base_pair_location < 39055519 + 1e6,]

#' Remove variants with very low frequency
hb = hb[minor_AF >=0.01 & minor_AF <=0.99,]
PCOS = PCOS[effect_allele_frequency >=0.01 & effect_allele_frequency <=0.99,]

#' format the information with rs ids and snps
head(hb)


#' Check overlap of rsID
hb[,rsID := gsub(":.*","",snp)]
PCOS = PCOS[rs_id %in% hb$rsID]
hb = hb[rsID %in% PCOS$rs_id, ]

# Check duplicates & triallelic SNPs
hb[duplicated(SNP),]
hb[duplicated(rsID),]
hb[duplicated(POS),]

#' Check order of data sets
setorder(PCOS,base_pair_location) #Alternative: match command - using id variable to order the other datasets 
stopifnot(hb$rsID == hb$rs_id)

#' # Harmonize alleles ####
#' ***
#EA and OA -> renaming the columns 
plot(hb$Freq_Tested_Allele, PCOS$effect_allele_frequency)

#' In the PCOS data, there are 3,383 SNPs with same allele coding as the BMI data, and 2,742 SNPs with switch alleles. These SNPs will be flipped back to the BMI coding. 
#' 
filt = hb$Tested_Allele == PCOS$other_allele & 
  hb$Other_Allele == PCOS$effect_allele
table(filt)
PCOS[filt,beta := beta * (-1)]
PCOS[filt,effect_allele_frequency := 1-effect_allele_frequency]
PCOS[filt,effect_allele := hb[filt,Tested_Allele]]
PCOS[filt,other_allele := hb[filt,Other_Allele]]

#' Now check the transformation with the EAF plot again. 
plot(BMI_comb$Freq_Tested_Allele, PCOS$effect_allele_frequency)

#' Remove the outlier that have different allele frquencies. 
filt = abs(BMI_fem$Freq_Tested_Allele - PCOS$effect_allele_frequency) >0.1
table(filt)
PCOS = PCOS[!filt,]
BMI_comb = BMI_comb[!filt,]


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

