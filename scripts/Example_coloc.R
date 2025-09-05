#' ---
#' title: "Colocalization"
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
#' Load data the harmonized data for BMI and check for shared signaling in males, females and the sex-combined data. 
#' 
#' All necessary R packages and paths to data files are in the source file
#' 
#' This is an classic R file that can be compiled into a report just like an Rmarkdown file. 
#' 
#' # Initialize ####
#' ***

rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load data ####
#' ***
load("../data/Input_harmonized.RData")

#' # Prepare data ####
#' ***
#' See Chris Wallace coloc vignette for more information: https://chr1swallace.github.io/coloc/articles/a02_data.html
#' 
#' ## Add missing columns 
#' 
#' We need MAF in all data sets and sample size in PCOS. 
#' 
BMI_fem[,MAF := ifelse(Freq_Tested_Allele<0.5,Freq_Tested_Allele,1-Freq_Tested_Allele)]
BMI_mal[,MAF := ifelse(Freq_Tested_Allele<0.5,Freq_Tested_Allele,1-Freq_Tested_Allele)]
BMI_comb[,MAF := ifelse(Freq_Tested_Allele<0.5,Freq_Tested_Allele,1-Freq_Tested_Allele)]
PCOS[,MAF := ifelse(effect_allele_frequency<0.5,effect_allele_frequency,1-effect_allele_frequency)]

PCOS[,Ncases := 14467]
PCOS[,Ncontrols := 430267]
PCOS[,Ntotal := 444734]

#' ## Create list objects 
#' 
#' Note: for PCOS, I use the position information of the BMI data, as PCOS is in hg38, not hg19. 
#' 
data_females = list(beta = BMI_fem$BETA,
                    varbeta = BMI_fem$SE^2,
                    snp = BMI_fem$rsID,
                    position = BMI_fem$POS,
                    type = "quant",
                    N = BMI_fem$N,
                    MAF = BMI_fem$MAF)

data_males = list(beta = BMI_mal$BETA,
                  varbeta = BMI_mal$SE^2,
                  snp = BMI_mal$rsID,
                  position = BMI_mal$POS,
                  type = "quant",
                  N = BMI_mal$N,
                  MAF = BMI_mal$MAF)

data_combined = list(beta = BMI_comb$BETA,
                     varbeta = BMI_comb$SE^2,
                     snp = BMI_comb$rsID,
                     position = BMI_comb$POS,
                     type = "quant",
                     N = BMI_comb$N,
                     MAF = BMI_comb$MAF)

data_PCOS = list(beta = PCOS$beta,
                 varbeta = PCOS$standard_error^2,
                 snp = PCOS$rs_id,
                 position = BMI_comb$POS,
                 type = "cc",
                 N = PCOS$Ntotal,
                 MAF = PCOS$MAF)

#' # Visualize locus ####
#' ***
plot_dataset(data_females)
plot_dataset(data_males)
plot_dataset(data_combined)
plot_dataset(data_PCOS) 

#' # Colocalization test ####
#' ***
#' See also https://chr1swallace.github.io/coloc/articles/a03_enumeration.html
my.res1 <- coloc.abf(dataset1=data_females,
                     dataset2=data_males)
sensitivity(my.res1,rule="H4 > 0.5") 

#' **Result females vs. males**: no shared signal (though there might be a second signal!)
#' 
my.res2 <- coloc.abf(dataset1=data_females,
                     dataset2=data_combined)
sensitivity(my.res2,rule="H4 > 0.5") 

#' **Result females vs. sex-combined data**: no shared signal (though there might be a second signal!)
#' 
my.res3 <- coloc.abf(dataset1=data_males,
                     dataset2=data_combined)
sensitivity(my.res3,rule="H4 > 0.5") 

#' **Result males vs. sex-combined data**: no shared signal (though there might be a second signal!)
#' 
my.res4 <- coloc.abf(dataset1=data_males,
                     dataset2=data_PCOS)
sensitivity(my.res4,rule="H4 > 0.5") 

#' **Result males vs. PCOS**: no shared signal, most likely only signal for BMI 
#' 
my.res5 <- coloc.abf(dataset1=data_combined,
                     dataset2=data_PCOS)
sensitivity(my.res5,rule="H4 > 0.5") 

#' **Result sex-combined vs. PCOS**: no shared signal, most likely only signal for BMI 
#' 
my.res6 <- coloc.abf(dataset1=data_females,
                     dataset2=data_PCOS)
sensitivity(my.res6,rule="H4 > 0.5") 

#' **Result females vs. PCOS**: no shared signal, most likely only signal for BMI 
#' 
#' # Save ####
#' ***
#' Save colocalization results
#' 
res = rbind(my.res1$summary,my.res2$summary,my.res3$summary,
            my.res3$summary,my.res5$summary,my.res6$summary)
coloc = as.data.frame(res)
setDT(coloc)
coloc[,trait1 := c("BMI females","BMI females","BMI males","BMI females","BMI males","BMI sex-combined")]
coloc[,trait2 := c("BMI females","BMI males","BMI sex-combined",rep("PCOS",3))]

save(coloc,file = "/Users/harshikamohanraj/Downloads/Coloc.RData")



### HbA1c -----------
hb[,MAF := ifelse(Freq_Tested_Allele<0.5,Freq_Tested_Allele,1-Freq_Tested_Allele)]
PCOS[,MAF := ifelse(effect_allele_frequency<0.5,effect_allele_frequency,1-effect_allele_frequency)]

PCOS[,Ncases := 14467]
PCOS[,Ncontrols := 430267]
PCOS[,Ntotal := 444734]

#' ## Create list objects 
#' 
#' Note: for PCOS, I use the position information of the BMI data, as PCOS is in hg38, not hg19. 
#' 

data_hb = list(beta = hb$BETA,
                     varbeta = hb$SE^2,
                     snp = hb$rsID,
                     position = hb$POS,
                     type = "quant",
                     N = hb$N,
                     MAF = hb$MAF)

data_PCOS = list(beta = PCOS$beta,
                 varbeta = PCOS$standard_error^2,
                 snp = PCOS$rs_id,
                 position = BMI_comb$POS,
                 type = "cc",
                 N = PCOS$Ntotal,
                 MAF = PCOS$MAF)

#' # Visualize locus ####
#' ***
plot_dataset(data_hb)
plot_dataset(data_PCOS) 

#' # Colocalization test ####
#' ***
my.res1 <- coloc.abf(dataset1=data_hb,
                     dataset2=data_PCOS)
sensitivity(my.res1,rule="H4 > 0.5") 

#' **Result HbA1c vs. PCOS**: no shared signal (though there might be a second signal!)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

