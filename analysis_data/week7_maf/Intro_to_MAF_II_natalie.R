# Natalie Fortunato and Sankalp Mrutyunjaya

# install and load maftools, TCGAbiolinks, and ggplot2 
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

# Set Path
setwd('/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/analysis_data')

# Read Clinical Data csv
clinical <- read.csv("brca_clinical_data.csv")

# Initialize maf_object
maf_query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
#GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)

# Initialize clincal_rad and clincal_drug
clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")

clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

# 1) Chosen variable: regimen_number
# Remove empty values
clinical_drug_therapy_cleaned_mask <- ifelse(clinical_drug$therapy_types == 'Chemotherapy'
                                             | clinical_drug$therapy_types == 'Hormone Therapy', 
                                             T, F)
clinical_drug_therapy_cleaned <- clinical_drug[clinical_drug_therapy_cleaned_mask,]

# 2) CoOncoplot
chemo_mask <- ifelse(clinical_drug_reg_cleaned$therapy_types == 'Chemotherapy', T, F)
hormone_mask <- ifelse(clinical_drug_reg_cleaned$therapy_types == 'Hormone Therapy', T, F)

chemo_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[chemo_mask]
hormone_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[hormone_mask]

chemo_maf <- subsetMaf(maf = maf_object, tsb = chemo_patient_barcodes)
hormone_maf <- subsetMaf(maf = maf_object, tsb = hormone_patient_barcodes)

coOncoplot(m1 = chemo_maf, m2 = hormone_maf, m1Name = 'Chemotherapy', 
           m2Name = 'Hormone Therapy', borderCol = NA)

#################### Answer question ##############################

# 3) Contingency Table
# Chosen gene: TP53

# Subset maf to TP53
TP53_maf <- subsetMaf(maf = maf_object, genes = 'TP53')

# TP53 Barcodes
TP53_gene_barcodes <- TP53_maf@clinical.data$Tumor_Sample_Barcode
num_TP53_pos <- length(TP53_gene_barcodes)

num_chemo <- length(chemo_patient_barcodes)
num_hormone <- length(hormone_patient_barcodes)

chemo_TP53_pos <- length(intersect(chemo_patient_barcodes, TP53_gene_barcodes))
hormone_TP53_pos <- length(intersect(hormone_patient_barcodes, TP53_gene_barcodes))
chemo_TP53_neg <- length(TP53_gene_barcodes) - chemo_TP53_pos
hormone_TP53_neg <- length(TP53_gene_barcodes) - hormone_TP53_pos

contig <- matrix(c(chemo_TP53_pos, chemo_TP53_neg, hormone_TP53_pos, hormone_TP53_neg),
                 nrow = 2)

mosaicplot(contig)

# 4) Colollipop Plot
lollipopPlot2(m1 = chemo_maf, 
              m2 = hormone_maf, 
              m1_name = 'Chemotherapy', 
              m2_name = 'Hormone Therapy', 
              gene = "TP53")


# 5) KM Plot
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == 'Alive', T, F)

mafSurvival(maf = maf_object,
            genes = "TP53", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)

