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

rm(list=ls()) #Remove any existing objects in R 
#setwd("C:/Users/Amy/Documents/MR_practicals/MR-Course2") #Change to your downloaded files location
coursedata = read.csv("/Users/harshikamohanraj/Downloads/coursedata.csv") #Load data
ratio.all = read.csv("/Users/harshikamohanraj/Downloads/summarized_data.csv", row=1) # Load answer data for later
attach(coursedata) #Attach coursedata to the R search path 
bx=ratio.all["bx",]
by=ratio.all["by",]
bxse=ratio.all["bxse",]
byse=ratio.all["byse",]

#Two-stage least squares and inverse-variance weighted estimates 
beta.ivw = sum(bx*by*byse^-2)/sum(bx^2*byse^-2)
beta.ivw #IVW estimate 
se.ivw   = 1/sqrt(sum(bx^2*byse^-2))
se.ivw #standard error of the IVW estimate 
#estimate from the 2SLS method was 0.570785 with a standard error of 0.2291837, which is very similar to the estimates from the IVW approach.


#Inverse variance weighted formula 
#meta-analysis of the ratio estimates from the individual variants, using the first-order standard errors.
#Calculate the ratio estimates and first-order standard errors for four genetic variants 
beta.ratio.all =  t(by/bx) 
se.ratio.all = t(byse/bx)
#meta-analysis command - fixed-effect estimate and standard error using the metagen command is identical to the IVW estimate and standard error
metagen(beta.ratio.all, se.ratio.all)
metagen(beta.ratio.all, se.ratio.all)$TE.fixed
metagen(beta.ratio.all, se.ratio.all)$seTE.fixed
# IVW method can also be motivated as a ratio estimate using a weighted allele score as an instrumental variable
score <- g1*as.numeric(bx[1]) + g2*as.numeric(bx[2]) + g3*as.numeric(bx[3]) + g4*as.numeric(bx[4])
#Calculate the ratio estimate and its standard error (first-order) using this score as an instrumental variable
ivmodel.score = ivreg(y~x|score, x=TRUE)
summary(ivmodel.score)
#fit the weighted linear regression model and obtain the causal estimate
BY<-t(by) # rotates data to a column vector
BX<-t(bx)
BYSE<-t(byse)
BXSE<-t(bxse)
regression<- lm(BY~BX-1, weights = BYSE^-2)
summary(regression) 
summary(regression)$coef[1]   
summary(regression)$coef[1,2]/summary(regression)$sigma 
#standard error is divided by the sigma quantity in the final line of code 

plot(BX, BY, xlim=c(min(BX-2*BXSE, 0), max(BX+2*BXSE, 0)),
     ylim=c(min(BY-2*BYSE, 0), max(BY+2*BYSE, 0)))
for (j in 1:length(BX)) {
  lines(c(BX[j],BX[j]), c(BY[j]-1.96*BYSE[j], BY[j]+1.96*BYSE[j]))
  lines(c(BX[j]-1.96*BXSE[j],BX[j]+1.96*BXSE[j]), c(BY[j], BY[j]))
}
abline(h=0, lwd=1); abline(v=0, lwd=1)
# add IVW estimate
abline(a=0, b=sum(bx*by*byseˆ(-2)/sum(bxˆ2*byseˆ(-2),col="red")

#lines around each point ->lines represent the 95% confidence intervals for the genetic associations with the exposure and with the outcome
detach(coursedata)
rm(list=ls())

#packages types - two sample MR = package from MR-Base automates the selection of SNPs and datasets 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MendelianRandomization)
ratio.all <- as.matrix(read.csv("/Users/harshikamohanraj/Downloads/summarized_data.csv", row=1))# Load answer data for later
coursedata = read.csv("/Users/harshikamohanraj/Downloads/coursedata.csv") #Load data
attach(coursedata) #Attach coursedata to the R search path 
bx.all=ratio.all["bx",]
by.all=ratio.all["by",]
bxse.all=ratio.all["bxse",]
byse.all=ratio.all["byse",]

#define a MRInput object to ensure the data is correctly formatted 
MRObject = mr_input(bx = bx.all, bxse = bxse.all, by = by.all, byse = byse.all)
#inverse-variance weighted (IVW) estimate
mr_ivw(MRObject)

#comparing fixed effects and random effects 
diabetes_data = read.csv("/Users/harshikamohanraj/Downloads/diabetes_data.csv")
detach(coursedata)
attach(diabetes_data)
MRObject2 = mr_input(bx = beta.exposure, bxse = se.exposure, 
                     by = beta.outcome, byse = se.outcome,snps = SNP)
mr_ivw(MRObject2, model="fixed")


#Graphs for the Inverse-variance weighted method
# scatter plot
mr_plot(MRObject2)

# single snp plot
mr_forest(MRObject2, ordered=TRUE)

# leave one out plot
mr_loo(MRObject2)

# funnel plot
mr_funnel(MRObject2)
     
#MR Egger and median-based methods using MR 
# creates mr_input objects
MR_LDLObject = mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse)
MR_HDLObject = mr_input(bx = hdlc, bxse = hdlcse, by = chdlodds, byse = chdloddsse)
# fit some egger and median models
mr_egger(MR_LDLObject)
mr_median(MR_LDLObject, weighting="weighted")
mr_median(MR_LDLObject, weighting="simple")
# graph mr models
mr_allmethods(MR_LDLObject)
mr_allmethods(MR_LDLObject, method="main")


# fit some egger and median models
mr_egger(MR_HDLObject)
mr_median(MR_HDLObject, weighting="weighted")
mr_median(MR_HDLObject, weighting="simple")
# graph mr models
mr_allmethods(MR_HDLObject)
mr_allmethods(MR_HDLObject, method="main")

mr_plot(MR_LDLObject)
mr_plot(MR_LDLObject, orientate=TRUE, line="egger")
mr_plot(MR_LDLObject, interactive=FALSE)
mr_plot(MR_LDLObject, interactive=FALSE, labels=TRUE)
mr_plot(MR_HDLObject)
mr_plot(MR_HDLObject, orientate=TRUE, line="egger")
mr_plot(MR_HDLObject, interactive=FALSE)
mr_plot(MR_HDLObject, interactive=FALSE, labels=TRUE)
   

#End of course


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

