#Supplementary tables 

tab1 <- as.table(rbind(c("Body mass index","Pulit SL","European", "806,834",	"https://zenodo.org/records/1251813"), c("Haemoglobin A1c",	"Neale lab", "Across 6 continental ancestry groups in the UK Biobank","344,182", "http://www.nealelab.is/uk-biobank/")))
dimnames(tab1) <-list(exposures = c("BMI", "HbA1c"),source = c("Phenotype", "Study", "Ancestry", "N cases (/controls)","Link"))
tab1

tab2 <- as.table(rbind(c("polycystic ovary syndrome (PCOS)","Venkatesh","European", "14,467/430,267",	"https://www.ebi.ac.uk/gwas/studies/GCST90483500"), c("All-cause mortality", "Timmers PR", "Parental longevity ", "607344/404896", "https://www.ebi.ac.uk/gwas/publications/30642433"), c("CAD", "CARDIoGRAMplusC4D Consortium", "European, South Asian, and East Asian","60,801/123,504", "https://cardiogramplusc4d.org/data-downloads/")))
dimnames(tab2) <-list(outcomes = c("PCOS", "All-cause mortality","CAD"),source = c("Phenotype", "Study", "Ancestry", "N cases (/controls)","Link"))
tab2

tab3 <- as.table(rbind(c("rs17757975", "6", "38214150 rs17757975:T:C","0.8462","0.0155","0.0024","1.57e-10","718688","0.995585","rs17757975","combined"), c("rs4714290","6","40003502 rs4714290:T:C","0.7046","0.0112","0.0019","2.439e-09","715538","0.996559","rs4714290","combined"), c("rs116208210","6", "39920827 rs116208210:T:C", "0.9878","0.0674","0.0122","2.977e-08","262817","0.963629","rs116208210","females"), c("rs2254336","6", "39032835 rs2254336:A:T","0.4763","-0.015","0.0025","3.325e-09","325594","0.986084", "rs2254336","males")))
dimnames(tab3) <-list(exposure = c("BMI_snp1", "BMI_snp2", "BMI_snp3", "BMI_snp4"),  summary_statistic = c("Exposure SNP", "CHR", "POS", "Freq_Tested_Allele", "BETA", "SE", "P", "N", "INFO", "rsID","setting"))
tab3

tab4 <- as.table(rbind(c("6", "38246374", "T", "C",  "-0.0036", "0.0344", "0.8589", "0.916", "rs17757975"), c("6", "40035763", "T", "C", "-0.0223", "0.057", "0.7154", "0.6959", "rs4714290"), c("6", "39953088", "T", "C", "0.2259", "0.2338", "0.9873", "0.334", "rs116208210"), c("6", "39065059", "A", "T", "0.0364", "0.0225", "0.4775", "0.106", "rs2254336")))
dimnames(tab4) <-list(exposure = c("PCOS_snp1", "PCOS_snp2", "PCOS_snp3", "PCOS_snp4"),  summary_statistic = c("chromosome", "POS", "EA", "OA", "beta", "standard_error", "Freq_Tested_Allele", "p_value", "rs_id"))
tab4


head(select(cad_1, rs_number, BP, reference_allele, other_allele , eaf, beta, se, p_value))
head(select(cad_2, rs_number, BP, reference_allele, other_allele , eaf, beta, se, p_value))
head(select(cad_3, rs_number, BP, reference_allele, other_allele , eaf, beta, se, p_value))
head(select(cad_4, rs_number, BP, reference_allele, other_allele , eaf, beta, se, p_value))


tab5 <- as.table(rbind(c("6", "38214150", "T", "C","0.859674","0.005345","0.009322","0.566436"), c("6", "40003502","T","C","0.707305", "0.003027", "0.007056", "0.667966"), c("6", "39920827", "T", "C", "0.986937", "-0.02457","0.03115","0.43025"), c("6","39032835","A","T","0.480492","0.010777","0.006491","0.09689")))
dimnames(tab5) <-list(exposure = c("cad_1", "cad_2", "cad_3", "cad_4"),  summary_statistic = c("rs_number", "BP", "EA", "OA", "eaf", "beta", "se", "p_value"))
tab5


