# Natalie Fortunato
# QBIO 490 Midterm Project

# Research Question: Is there a correlation between the treatment type (hormonal or 
# chemotherapy) and the regulation of RYR2? 


# Set Working Directory to Midterm Project Folder
setwd("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/midterm_project_fortunato")

# Install and Load Libraries
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

# Query Data (query, download, prepare)
clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

maf_query <- GDCquery(
  project = "TCGA-BRCA", 
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

rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

# Create "outputs" Folder
dir.create('outputs')
setwd('outputs')

#############################
# Scatterplot of TP53 vs. RYR2 Gene Counts
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
rownames(rna_genes) <- rna_genes$gene_id

rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)

RYR2_mask <- ifelse(rna_genes$gene_name == "RYR2", T, F)
TP53_mask <- ifelse(rna_genes$gene_name == "TP53", T, F)

RYR2_counts <- rna_counts[RYR2_mask, ]
TP53_counts <- rna_counts[TP53_mask, ]


plot(as.numeric(RYR2_counts[1,]),
     as.numeric(TP53_counts[1,]),
     xlab = "RYR2 Counts",
     ylab = "TP53 Counts",
     main = "TP53 vs. RYR2 Counts",
     log = "xy")

# Save Plot
# Open File
jpeg("gene_counts_scatter.jpg")

# Create Plot
plot(as.numeric(RYR2_counts[1,]),
     as.numeric(TP53_counts[1,]),
     xlab = "RYR2 Counts",
     ylab = "TP53 Counts",
     main = "TP53 vs. RYR2 Counts",
     log = "xy")

# Close File
dev.off()

#############################
# KM Plot of Hormone vs. Chemotherapy

clinical_drug$Tumor_Sample_Barcode <- clinical_drug$bcr_patient_barcode
clinic_drug_merge <- merge(clinic, clinical_drug, by = "Tumor_Sample_Barcode")

# Remove Patients who did not Receive Hormone or Chemotherapy
therapy_mask <- ifelse(clinic_drug_merge$therapy_types == 'Chemotherapy'
                          | clinic_drug_merge$therapy_types == 'Hormone Therapy',
                          T, F)
therapy_cleaned_clinic <- clinic_drug_merge[therapy_mask,]

# Remove Patients without vital_status Data
na_mask <- ifelse(therapy_cleaned_clinic$vital_status == "<NA>", F, T)
na_mask <- !is.na(na_mask)
therapy_cleaned_clinic <- therapy_cleaned_clinic[na_mask,]

# Make a Survival Time Column for Survival Plots
therapy_cleaned_clinic$survival_time <- ifelse(is.na(therapy_cleaned_clinic$days_to_death),
                                           therapy_cleaned_clinic$survival_time
                                           <- therapy_cleaned_clinic$days_to_last_followup,
                                           therapy_cleaned_clinic$survival_time 
                                           <- therapy_cleaned_clinic$days_to_death)

# Remove Any Infinite Values in Survival Time Column
inf_mask <- ifelse(therapy_cleaned_clinic$survival_time == "-Inf", F, T)
therapy_cleaned_clinic <- 
  therapy_cleaned_clinic[inf_mask, ]

# Make a Death Event (T/F) Column for Survival Plots
therapy_cleaned_clinic$death_event <- ifelse(therapy_cleaned_clinic$vital_status == "Alive", F, T)

# Initialize a Survival Object
surv_object_therapy <- Surv(time = therapy_cleaned_clinic$survival_time, event = 
                          therapy_cleaned_clinic$death_event)

# Create a Fit Object 
therapy_fit <- survfit(surv_object_therapy ~ therapy_cleaned_clinic$therapy_types, data = therapy_cleaned_clinic)

# Format and Create KM plot 
survplot_therapy <- ggsurvplot(therapy_fit, pval=TRUE, 
                               ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                               legend = "right")

# Save KM Plot to a Variable
KM_plot_therapy <- survplot_therapy$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))

# Show Plot
KM_plot_therapy

# Save Plot
jpeg("KM_plot_therapy.jpg")
KM_plot_therapy <- survplot_therapy$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))
KM_plot_therapy
dev.off()


#############################
# Lollipop Plot of RYR2 Mutations in Hormone vs. Chemotherapy

# Create a Vector of Chemo/Horomone Patient Barcodes and Use subsetMaf to create
# a MAF with only Chemo/Hormone Patients
chemo_mask <- ifelse(clinic_drug_merge$therapy_types == 'Chemotherapy', T, F)
chemo_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[chemo_mask]
chemo_maf <- subsetMaf(maf = maf_object,
                       tsb = chemo_patient_barcodes)

hormone_mask <- ifelse(clinic_drug_merge$therapy_types == 'Hormone Therapy', T, F)
hormone_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[hormone_mask]
hormone_maf <- subsetMaf(maf = maf_object,
                       tsb = hormone_patient_barcodes)

# Save Plot to Variable
lollipop_therapy_plot <- lollipopPlot2(m1 = chemo_maf, 
                                      m2 = hormone_maf, 
                                      m1_name = 'Pateints who Recieved Chemotherapy',
                                      m2_name = 'Patients who Recieved Hormone Therapy',
                                      gene = "RYR2")

# Show Plot
lollipop_therapy_plot

# Save Plot
jpeg("lollipop_therapy_plot.jpg")
lollipop_therapy_plot <- lollipopPlot2(m1 = chemo_maf, 
                                       m2 = hormone_maf, 
                                       m1_name = 'Pateints who Recieved Chemotherapy',
                                       m2_name = 'Patients who Recieved Hormone Therapy',
                                       gene = "RYR2")

lollipop_therapy_plot
dev.off()


#############################
# Volcano Plot

# Load DESeq2 Library
library(DESeq2)

# Data Preparation
# Create rna_clinical which Contains Data from colData
rna_clinical <- rna_se@colData[!is.na(rna_se@colData$age_at_index), ]
rna_clinical <- as.data.frame(rna_clinical)

# Subset out Treatments Column
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F, T)
rna_clinical <- rna_clinical[ ,treatments_mask]

rna_counts_age <- rna_se@assays@data$unstranded[, !is.na(rna_se@colData$age_at_index)]
rna_counts_age <- as.data.frame(rna_counts_age)

# Name rna_counts_age Rows with gene_id
rownames(rna_counts_age) <- rna_genes$gene_id

# Create Column with Age Factor
rna_clinical$age_category <- factor(ifelse(rna_clinical$age_at_index <= 58, "young", "old"))

# Mask Out "Solid Tissue Normal"
tumor_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_counts_age <- rna_counts_age[ ,tumor_mask]
rna_clinical <- rna_clinical[tumor_mask, ]
rna_clinical$definition <- factor(rna_clinical$definition)

# Remove Genes where total gene_counts across all Patients < 10
gene_counts_mask <- ifelse(rowSums(rna_counts_age) < 10, F, T)
rna_counts_age <- rna_counts_age[gene_counts_mask, ]

# Make Sure colnames of rna_counts_age and rownames of rna_clinical Match
colnames(rna_counts_age) <- rownames(rna_clinical)

# Create DESeq Data Set dds
dds <- DESeqDataSetFromMatrix(countData = rna_counts_age,
                              colData = rna_clinical,
                              design = ~ definition)

# Create DESeq Object
dds_obj <- DESeq(dds)

# See Results --------
resultsNames(dds_obj)

# See Young vs. Old Comparison
results <- results(dds_obj, format = "DataFrame",
                   contrast = c("definition", "Primary solid Tumor", "Metastatic"))

gene_id_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
rna_genes <- rna_genes[gene_id_mask,]

# Save Results to a Variable
results <- data.frame(rna_genes$gene_name, rownames(results),
                      results$log2FoldChange, results$pvalue, results$padj,
                      -log10(results$padj))

# Load EnhancedVolcano Library
library(EnhancedVolcano)

# Save Plot to a Variable
volcano_plot <-   EnhancedVolcano(results,
                                  lab = rownames(results),
                                  x = 'results.log2FoldChange',
                                  y = 'results.pvalue')

# Show Plot
volcano_plot

# Save Plot
jpeg("volcano_plot.jpg")
volcano_plot <-   EnhancedVolcano(results,
                                  lab = rownames(results),
                                  x = 'results.log2FoldChange',
                                  y = 'results.pvalue')
volcano_plot
dev.off()
