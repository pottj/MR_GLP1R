#' ---
#' title: "Example script for colocalization"
#' subtitle: ""
#' author: "Janne Pott"
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
#' In this example script, I load data the harmonized data for BMI and check for shared signaling in males, females and the sex-combined data. 
#' 
#' All necessary R packages and paths to data files are in the source file
#' 
#' This is an classic R file that can be compiled into a report just like an Rmarkdown file. 
#' 
#' # Initialize ####
#' ***

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ivreg, meta, MendelianRandomization)

rm(list=ls()) #Remove any existing objects in R 
pacman::p_load("car")
library(car)
attach(Prestige) #Attach Prestige to the R search path 
help(Prestige) #Description of the dataset 
str(Prestige) #Info on the structure of the data 
head(Prestige) #Show first 6 entries 
plot(education, income)
model<-lm(income ~ education)
summary(model)
edu_income = summary(model) #Store results from the regression model 
slope = edu_income$coef[2] #Estimate of the slope parameter  
slope 
slope_se = edu_income$coef[2,2] #Standard error of the slope parameter  
slope_se
f_stat = edu_income$f[1] #F-statistic 
f_stat 
PredictedValues<- predict(model) # the fitted values for the model data
cbind(Actual = income[1:5], Fitted= model$fitted.values[1:5],  Predicted= PredictedValues[1:5])
new <- data.frame(education = seq(5, 15, 2)) #predicted values for a new data set
PredictedValues2<- predict(model, newdata=new, se.fit = TRUE)
cbind(Education=new, PredictedIncome=PredictedValues2$fit, SE= PredictedValues2$se.fit)
plot(education, income)
abline(model, col="red")
summary(lm(income ~ education-1))


rm(list=ls()) #Remove any existing objects in R 
#setwd("C:/Users/Amy/Documents/MR_practicals/MR-Course2") #Change to your downloaded files location
coursedata = read.csv("/Users/harshikamohanraj/Downloads/coursedata.csv") #Load data
ratio.all = read.csv("/Users/harshikamohanraj/Downloads/summarized_data.csv", row=1) # Load answer data for later
attach(coursedata) #Attach coursedata to the R search path 
str(coursedata) #Info on the structure of the data 
head(coursedata) #Show first 6 entries 
#casual estimate using the ratio method for a continuous outcome 
#x=risk factor, y=continuous outcome 
#calculate the ratio causal estimate for the first genetic variant g1
by1 = lm(y~g1)$coef[2] #Genetic association with the outcome
bx1 = lm(x~g1)$coef[2] #Genetic association with the exposure
beta.ratio1 = by1/bx1  
beta.ratio1 #Ratio estimate for g1 
#The standard error of the causal estimate = SE of genetic association with the outcome/genetic association with the risk factor 
byse1 = summary(lm(y~g1))$coef[2,2] #Standard error of the G-Y association
se.ratio1first = byse1/sqrt(bx1^2)  
se.ratio1first #Standard error (first order) of the ratio estimate 
#above approximation does not account for the uncertainty in the denominator of the ratio estimate
bxse1 = summary(lm(x~g1))$coef[2,2] #Standard error of the G-X association
se.ratio1second = sqrt(byse1^2/bx1^2 + by1^2*bxse1^2/bx1^4)
se.ratio1second #Standard error (second order) of the ratio estimate 
# F-statistic from the regression of the risk factor on the genetic variant(s) is used as a measure of ‘weak instrument bias’
fstat1 = summary(lm(x~g1))$f[1]
fstat1 #Some studies recommend excluding genetic variants if they have a F-statistic less than 10
#Minor Allele Frequency (MAF) is the frequency at which the second most common allele occurs in a given population
#MAF for g1
MAF = (sum(g1==1) + 2*sum(g1==2))/(2*length(g1))
MAF
ratio.all

#observational association by regression y on x 
lm(y~x)$coef[2]

#two-stage least squares method “by hand”
by.hand.fitted.values<-lm(x~g1+g2+g3+g4)$fitted
by.hand<-lm(y~by.hand.fitted.values)
summary(by.hand)$coef[2]  #estimate
summary(by.hand)$coef[2,2] # standard error

#Ivreg function for two stage least squares method 
ivmodel.all = ivreg(y~x|g1+g2+g3+g4, x=TRUE)
summary(ivmodel.all)$coef[2] #2SLS estimate 
summary(ivmodel.all)$coef[2,2] #Standard error of the 2SLS estimate 
#The estimates are the same, but the standard error is slightly larger using the ivreg function as it takes into account the uncertainty in the first stage regression.

#F statistic
summary(lm(x~g1+g2+g3+g4))$f[1]
#the two-stage least squares method based only on the first genetic variant g1.
ivmodel.g1 = ivreg(y~x|g1, x=TRUE)
summary(ivmodel.g1)$coef[2] #2SLS estimate for g1 only 
summary(ivmodel.g1)$coef[2,2] #Standard error of the 2SLS estimate for g1 only

#Causal estimate for a binary outcome 
#causal effect of the risk factor x on the binary outcome y.bin -> estimated as the per allel genetic association with the outcome divided  by the per allele genetic association with the risk factor 
#with a binary outcome, logistic regression used for regression of the outcome on the genetic variant 
#Evaluate the ratio estimate and its standard error for g1 using logistic regression
by1.bin   = glm(y.bin~g1, family=binomial)$coef[2] #logistic regression for - numerator (gene-outcome association)
#gene-outcome association
byse1.bin = summary(glm(y.bin~g1, family=binomial))$coef[2,2]
bx1.bin   = lm(x[y.bin==0]~g1[y.bin ==0])$coef[2] #linear regression in the controls only - denominator calculate 
beta.ratio1.bin = by1.bin/bx1.bin
beta.ratio1.bin #ratio estimate for g1
se.ratio1.bin   = byse1.bin/bx1.bin
se.ratio1.bin #standard error of the ratio estimate for g1

#Calculate the two-stage least squares estimate for g1 only
#first regress the risk factor on the genetic variant - regression only on the controls 
#then, find fitted values for the model 
g1.con = g1[y.bin ==0] #values for g1 in the controls only 
x.con  = x[y.bin ==0] #values for the risk factor in the controls only 
predict.con.g1 = predict(lm(x.con~g1.con), newdata=list(g1.con=g1)) #Generate predicted
#values for all participants based on the linear regression in the controls only.  
tsls1.con = glm(y.bin~predict.con.g1, family=binomial) #Fit a logistic regression 
#model on all the participants
summary(tsls1.con)$coef[2]
summary(tsls1.con)$coef[2,2] 
#two-stage least squares estimate for all of the genetic variants.
g2.con = g2[y.bin ==0] #values for g2 in the controls only 
g3.con = g3[y.bin ==0] #values for g3 in the controls only 
g4.con = g4[y.bin ==0] #values for g4 in the controls only
#Predicted values 
predict.con<-predict(lm(x.con~g1.con+g2.con+g3.con+g4.con), newdata=c(list(g1.con=g1),list(g2.con=g2), list(g3.con=g3),list(g4.con=g4)))
tsls1.con.all = glm(y.bin~predict.con, family=binomial) #Logistic regression
summary(tsls1.con.all)$coef[2] #
summary(tsls1.con.all)$coef[2,2] # represents the change in log causal odds ratio for y.bin per one unit increase in x.

#cleanup workspace 
detach(coursedata)
rm(list=ls())


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
                 snp = PCOS$rsid,
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

#' **Result males vs. sex-combined data**: shared signal (though there might be a second signal!)
#' 
my.res4 <- coloc.abf(dataset1=data_females,
                     dataset2=data_PCOS)
sensitivity(my.res4,rule="H4 > 0.5") 

#' **Result males vs. PCOS**: no shared signal, most likely only signal for BMI 
#' 
my.res5 <- coloc.abf(dataset1=data_males,
                     dataset2=data_PCOS)
sensitivity(my.res5,rule="H4 > 0.5") 

#' **Result sex-combined vs. PCOS**: no shared signal, most likely only signal for BMI 
#' 
my.res6 <- coloc.abf(dataset1=data_combined,
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

save(coloc,file = "../results/Coloc.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

