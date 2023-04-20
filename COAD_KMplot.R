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
#clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
#clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

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
#colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

# Determine if there are any NA values (there are none)
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))

#is_na_age <- is.na(clinic)

# Create a column in clinic with Young/Old 
young_mask <- ifelse((clinic$age_at_initial_pathologic_diagnosis < 50), T, F)
clinic$age_status <- ifelse(young_mask, "Young", "Old")


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

#initializing a survival object
surv_object_age <- Surv(time = cleaned_clinic$survival_time, event = 
                          cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ cleaned_clinic$age_status, data = cleaned_clinic)

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







