# Natalie Fortunato

setwd("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/analysis_data")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")

library(BiocManager)

library(TCGAbiolinks)

library(maftools)

clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

sum(is.na(clinic$age_at_initial_pathologic_diagnosis))
# 1) I chose the variable age_at_initial_pathologic_diagnosis. It is the age that the patient
# was when they were diagnosed.
# 2) This variable is discrete.

sum(is.na(clinical_rad$radiation_dosage))
# 3) I chose the variable radiation_dosage. 
# 4) This variable is discrete.

# 5) a. Patients diagnosed at the same age were treated with similar radiation dosages.
#    b. Patients diagnosed at a younger age have a higher rate of survival.
#    c. Pateints treated with a higher dosage of radiation have a higher rate of survival.

clinic_rad_merge <- merge(clinic, clinical_rad, by = "bcr_patient_barcode")

# create scatterplot
plot(x = clinic_rad_merge$age_at_initial_pathologic_diagnosis, 
     y = clinic_rad_merge$radiation_dosage, 
     main = 'Radiation Dosage vs. Age at Diagnosis', 
     xlab = 'Age at Initial Pathologic Diagnosis', ylab = 'Radiation Dosage')

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/analysis_data/
     scatterplot_rad_dose_age.jpg")
scatterplot_rad_dose_age <- plot(x = clinic_rad_merge$age_at_initial_pathologic_diagnosis, 
                                 y = clinic_rad_merge$radiation_dosage, 
                                 main = 'Radiation Dosage vs. Age at Diagnosis', 
                                 xlab = 'Age at Initial Pathologic Diagnosis', 
                                 ylab = 'Radiation Dosage')
scatterplot_rad_dose_age
dev.off()

# 1) From the plot I can see that 80 and 100 were very common radiation dosages, but there
#    doesn't seem to be a correlation between the patient's radiation dosage and their
#    age at the time of diagnosis. I chose to create a scatterplot to demonstrate the 
#    relationship between these two variables because they are both discrete variables and
#    not categorical.

if (!require(survival)){
  install.packages("survival")
}
library(survival)

if (!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if (!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)


####### Age KM Plot
# remove patients with no age info
age_na_mask <- ifelse(is.na(clinic_rad_merge$age_at_initial_pathologic_diagnosis), F, T)
age_cleaned_clinic <- clinic[age_na_mask, ]

# create a column in age_cleaned_clinic where patients are labeled "Young", "Middle", or "Old
young_mask <- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis < 35, T, F)
middle_mask <- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis >= 35 
                      & age_cleaned_clinic$age_at_initial_pathologic_diagnosis <= 60, T, F)
old_mask <- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis >60, T, F)
age_cleaned_clinic$age_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle",
                                        "Old"))

#making a survival time column for survival plots
age_cleaned_clinic$survival_time <- ifelse(is.na(age_cleaned_clinic$days_to_death),
                                           age_cleaned_clinic$survival_time
                                           <- age_cleaned_clinic$days_to_last_followup,
                                           age_cleaned_clinic$survival_time 
                                           <- age_cleaned_clinic$days_to_death)

#remove any -Inf values in survival_time
inf_mask <- ifelse(age_cleaned_clinic$survival_time == "-Inf", F, T)
age_cleaned_clinic <- age_cleaned_clinic[inf_mask, ]

#making a death event (T/F) column for survival plots
age_cleaned_clinic$death_event < ifelse(age_cleaned_clinic$vital_status == "Alive",
                                        age_cleaned_clinic$death_event <- FALSE,
                                        age_cleaned_clinic$death_event <- TRUE)

#initializing a survival object
surv_object_age <- Surv(time = age_cleaned_clinic$survival_time, event = 
                          age_cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ age_cleaned_clinic$age_status, data = age_cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))

#show plot
KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/analysis_data/KM_plot_age.jpg")
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))
KM_plot_age
dev.off()


############ Radiation Dosage KM Plot
# remove patients with no radiation info
rad_na_mask <- ifelse(is.na(clinical_rad$radiation_dosage), F, T)
cleaned_clinical_rad <- clinical_rad[rad_na_mask, ]

#change radiation_dosage from factor to numeric type
cleaned_clinical_rad$radiation_dosage <- as.numeric(cleaned_clinical_rad$radiation_dosage)

# create a column in age_cleaned_clinic where patients are labeled "Low", "Medium", or "High"
low_mask <- ifelse(cleaned_clinical_rad$radiation_dosage < 40, T, F)
medium_mask <- ifelse(cleaned_clinical_rad$radiation_dosage >= 40 
                      & cleaned_clinical_rad$radiation_dosage <= 80, T, F)
high_mask <- ifelse(cleaned_clinical_rad$radiation_dosage > 80, T, F)
cleaned_clinical_rad$rad_status <- ifelse(low_mask, "Low", ifelse(medium_mask, "Medium",
                                                                    "High"))


#clinic_merge <- merge(clinic_rad, clinical, by = "bcr_patient_barcode")

#making a survival time column for survival plots
clinic_rad_merge$survival_time <- ifelse(is.na(clinic_rad_merge$days_to_death),
                                             clinic_rad_merge$survival_time
                                           <- clinic_rad_merge$days_to_last_followup,
                                           clinic_rad_merge$survival_time 
                                           <- clinic_rad_merge$days_to_death)

#remove any -Inf values in survival_time
inf_mask <- ifelse(clinic_rad_merge$survival_time == "-Inf", F, T)
clinic_rad_merge <- clinic_rad_merge[inf_mask, ]

#making a death event (T/F) column for survival plots
clinic_rad_merge$death_event < ifelse(clinic_rad_merge$vital_status == "Alive",
                                      clinic_rad_merge$death_event <- FALSE,
                                      clinic_rad_merge$death_event <- TRUE)

#initializing a survival object
surv_object_rad <- Surv(time = clinic_rad_merge$survival_time, event = 
                          clinic_rad_merge$death_event)

#creating a fit object
rad_fit <- survfit(surv_object_rad ~ clinic_rad_merge$vital_status, data = clinic_rad_merge)

#format and create KM plot
survplot_rad <- ggsurvplot(rad_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
KM_plot_rad <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))

#show plot
KM_plot_rad

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/analysis_data/KM_plot_rad.jpg")
KM_plot_rad <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))
KM_plot_rad
dev.off()

# 4) Analyze your two KM plots. What do the KM plots suggest about the impact of 
#    your variables on breast cancer survival? What are the p-values? Do the 
#    differences in survival appear to be significant?

# The age KM plot had a p value of <0.0001 which means that the results were
# significant suggesting that patients diagnosed when they were young had
# a higher probability of survival, followed by middle and lastly old.
# The 

