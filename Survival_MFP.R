library(dplyr)
library(ggplot2)
library(gtsummary)
library(kableExtra)
library(survival)
library(survminer)

mfp = read.table("../TCGA_MFP/MFP_labels_TCGA.tsv",h=T,sep="\t")
colnames(mfp)[1] = "PatientID"

#clinical
clinHNSC <- read.table("hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",h=T,sep="\t")
clinHNSC <- clinHNSC[,c("Patient.ID","Sample.ID","Overall.Survival.Status","Overall.Survival..Months.",
                    "Disease.Free..Months.","Progress.Free.Survival..Months.")]
colnames(clinHNSC)[1] = "PatientID"
allHNSC <- merge(clinHNSC,mfp,by='PatientID')

#particular case TNBC
clinTNBC <- read.table("brca_tcga_pan_can_atlas_2018_clinical_data.tsv",h=T,sep="\t")
clinTNBC <- clinTNBC[,c("Patient.ID","PatientID_other","Sample.ID","Overall.Survival.Status","Overall.Survival..Months.",
                        "Disease.Free..Months.","Progress.Free.Survival..Months.")]
tnbc <- read.table("TNBC_cases.txt",h=T,sep="\t")
clinTNBC <- clinTNBC[clinTNBC$Patient.ID %in% tnbc$PatientID,]
colnames(clinTNBC)[1] = "PatientID"
allTNBC <- merge(clinTNBC,mfp,by='PatientID')

#Survival curves
fitMFP <- survfit(Surv(Overall.Survival..Months., 
                       Overall.Survival.Status) ~MFP, 
                  data = allHNSC)
ggsurvplot(fitMFP, 
           data = allHNSC, 
           risk.table = TRUE,
           palette = c("#bd8c2b","#bc2a46","#493ac2","#6ac0c1"),
           pval = TRUE,
           font.x = 16,
           font.y = 16,
           title="Overall Survival using MFP classification",
           conf.int = FALSE)

fitPFS <- survfit(Surv(Progress.Free.Survival..Months., 
                       Overall.Survival.Status) ~MFP, 
                  data = allHNSC)
ggsurvplot(fitPFS, 
           data = allHNSC, 
           risk.table = TRUE,
           palette = c("#bd8c2b","#bc2a46","#493ac2","#6ac0c1"),
           pval = TRUE,
           font.x = 16,
           font.y = 16,
           title="Progress-Free Survival using MFP classification",
           conf.int = FALSE)
