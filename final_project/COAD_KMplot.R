# Natalie Fortunato
# Colorectal Cancer Kaplan Meier Plot

# Set Working Directory
setwd("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato")


# Download Libraries
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
library(survival)
library(survminer)
library(ggplot2)

# Query Data (query, download, prepare)
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")

# Error in TCGABiolinks Query (missing vital_status data)
clinic <- read.csv("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/harmonized_coad_clinical.csv")

maf_query <- GDCquery(
  project = "TCGA-COAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

# Change Column Names of Clinic to Patient Barcodes
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

# Determine if there are any NA values (there are none)
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))


# Create a column in clinic with Young/Old 
young_mask <- ifelse((clinic$age_at_initial_pathologic_diagnosis < 50), T, F)
clinic$age_status <- ifelse(young_mask, "Early-Onset", "Late-Onset")


#making a survival time column for survival plots
clinic$survival_time <- ifelse(is.na(clinic$days_to_death),
                               clinic$survival_time <- clinic$days_to_last_followup,
                               clinic$survival_time <- clinic$days_to_death)

#remove any -Inf and NA values in survival_time
cleaned_clinic <- clinic
na_mask <- ifelse(is.na(clinic$survival_time), F, T)
cleaned_clinic <- cleaned_clinic[na_mask, ]
inf_mask <- ifelse(cleaned_clinic$survival_time == "-Inf", F, T)
cleaned_clinic <- cleaned_clinic[inf_mask, ]

#making a death event (T/F) column for survival plots
cleaned_clinic$death_event <- ifelse(cleaned_clinic$vital_status == "Alive",
                                    cleaned_clinic$death_event <- FALSE,
                                    cleaned_clinic$death_event <- TRUE)

################ Male Plot ####################
male_mask <- ifelse((cleaned_clinic$gender == 'MALE'), T, F)
male_cleaned_clinic <- cleaned_clinic[male_mask, ]

#initializing a survival object
surv_object_age <- Surv(time = male_cleaned_clinic$survival_time, event = 
                          male_cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ male_cleaned_clinic$age_status, data = male_cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
male_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))

#show plot
male_KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/male_COAD_KM_plot_age.jpg")
male_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))
male_KM_plot_age
dev.off()


################ Female Plot ####################
female_mask <- ifelse((cleaned_clinic$gender == 'FEMALE'), T, F)
female_cleaned_clinic <- cleaned_clinic[female_mask, ]

#initializing a survival object
surv_object_age <- Surv(time = female_cleaned_clinic$survival_time, event = 
                          female_cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ female_cleaned_clinic$age_status, data = female_cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
female_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                           axis.text = element_text(size=16),
                                                           legend.title = element_text(size=14),
                                                           legend.text = element_text(size=12))

#show plot
female_KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/female_COAD_KM_plot_age.jpg")
female_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                           axis.text = element_text(size=16),
                                                           legend.title = element_text(size=14),
                                                           legend.text = element_text(size=12))
female_KM_plot_age
dev.off()






