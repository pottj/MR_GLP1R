#Supplementary tables 

tab1 <- as.table(rbind(c("Body mass index","Pulit SL","European", "806,834",	"https://zenodo.org/records/1251813"), c("Haemoglobin A1c",	"Neale lab", "Across 6 continental ancestry groups in the UK Biobank","344,182", "http://www.nealelab.is/uk-biobank/"), c("Coronary artery disease", "")))
dimnames(tab) <-list(exposure = c("BMI", "HbA1c"),source = c("Phenotype", "Study", "Ancestry", "N cases (/controls)","Link"))
tab

tab2 <- as.table(rbind(c("polycystic ovary syndrome (PCOS)","Venkatesh","European", "14,467/430,267",	"https://www.ebi.ac.uk/gwas/studies/GCST90483500"), c("All-cause mortality", "Timmers PR", "Parental longevity ", "607344/404896", "https://www.ebi.ac.uk/gwas/publications/30642433"), c("CAD", "CARDIoGRAMplusC4D Consortium", "European, South Asian, and East Asian","60,801/123,504", "https://cardiogramplusc4d.org/data-downloads/")))
dimnames(tab2) <-list(exposure = c("PCOS", "All-cause mortality","CAD"),source = c("Phenotype", "Study", "Ancestry", "N cases (/controls)","Link"))
tab2

