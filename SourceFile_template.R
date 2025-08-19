#############################
# this is the template source file 
#############################

#############################
# R library and R packages
#############################
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MendelianRandomization))
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(WriteXLS))

#############################
# path to data (root is the scripts folder)
#############################
# Please update the path to your data here!

# GWAS summary statistics for BMI

Pulit_BMI_female = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit_BMI_male = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit_BMI_combined = "/Users/harshikamohanraj/Downloads/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt"

# GWAS summary statistics for PCOS

Venkatesh_PCOS = "Users/harshikamohanraj/Downloads/CST90483500.h.tsv.gz"

#############################
# Helper function
#############################
MRfunction_jp<-function(betaX,seX,betaY,seY){
  betaIV<-betaY/betaX
  se.ratio.st<- seY/sqrt(betaX^2)
  se.ratio.nd<- sqrt(seY^2/betaX^2 + betaY^2*seX^2/betaX^4)
  p1<-2*pnorm(-abs(betaIV/se.ratio.st))
  p2<-2*pnorm(-abs(betaIV/se.ratio.nd))
  FStat <- (betaX/seX)^2
  res<-data.table(
    beta_IV=betaIV,
    se_IV1=se.ratio.st,
    se_IV2=se.ratio.nd,
    p_IV1=p1,
    p_IV2=p2, 
    Fstat=FStat)
  return(res)
}
