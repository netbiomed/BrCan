#================================================================
#
#
# Author: Mauricio Reyes-Elizondo
# Last modified: Jul 03, 2024
# Project: BrCan
# Place: Instituto Nacional de Medicina Genómica (CDMX, Remote)
#
#
# ===============================================================

# __        __         _    _                
# \ \      / /__  _ __| | _(_)_ __   __ _    
#  \ \ /\ / / _ \| '__| |/ / | '_ \ / _` |   
#   \ V  V / (_) | |  |   <| | | | | (_| |   
#    \_/\_/ \___/|_|  |_|\_\_|_| |_|\__, |   
#   __| (_)_ __ ___  ___| |_ ___  _ |___/  _ 
#  / _` | | '__/ _ \/ __| __/ _ \| '__| | | |
# | (_| | | | |  __/ (__| || (_) | |  | |_| |
#  \__,_|_|_|  \___|\___|\__\___/|_|   \__, |
#                                      |___/ 

setwd("~/Documents/BrCanProject/") # Change directory

load("~/Documents/BrCanProject/metadata_image.RData") # Change directory

# ===============================================================

#  _     _ _                    _           
# | |   (_) |__  _ __ __ _ _ __(_) ___  ___ 
# | |   | | '_ \| '__/ _` | '__| |/ _ \/ __|
# | |___| | |_) | | | (_| | |  | |  __/\__ \
# |_____|_|_.__/|_|  \__,_|_|  |_|\___||___/

library(TCGAbiolinks)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)
library(maftools)

# ===============================================================

#  ____      _        _                  ____       _   _            _       
# |  _ \ ___| |_ _ __(_) _____   _____  |  _ \ __ _| |_(_) ___ _ __ | |_ ___ 
# | |_) / _ \ __| '__| |/ _ \ \ / / _ \ | |_) / _` | __| |/ _ \ '_ \| __/ __|
# |  _ <  __/ |_| |  | |  __/\ V /  __/ |  __/ (_| | |_| |  __/ | | | |_\__ \
# |_| \_\___|\__|_|  |_|\___| \_/ \___| |_|   \__,_|\__|_|\___|_| |_|\__|___/


# Retrieve TCGA-BRCA patients with transcriptomic data
brca <- curatedTCGAData(
  diseaseCode = "BRCA", assays = c("Mutation", "RNASeq2*"), 
  version = "1.1.38", 
  dry.run = FALSE
)

sampleTables(brca)
# Updated as of July 3, 2024
#
# $`BRCA_Mutation-20160128`
#
# 01  06 
# 988   5 
#
# $`BRCA_RNASeq2GeneNorm-20160128`
#
# 01   06   11 
# 1093    7  112 

brca_samples <- as.data.frame(brca@sampleMap) # Convert data as a data frame
# write.csv(brca_samples, "brcacodes_moleculardata.csv", row.names = TRUE)
# saveRDS(brca_samples, "brcacodes_moleculardata.RDS")
head(brca_samples)

brca_samples$sample_types <- substring(brca_samples$colname, 14, 15) # Divide samples
unique(brca_samples$sample_types) # "01" "01" "06"

# 01: Primary tumor
# 06: Metastastic tumor
# 11: Adjacent normal tissue

brca_samples <- brca_samples[!(brca_samples$sample_types == "06"),] # Remove metastasis samples

# Retrieve patients with available transcriptomic data
brca_rna_codes <- brca_samples %>% 
  filter(assay == "BRCA_RNASeq2GeneNorm-20160128", sample_types == "01") %>%
  pull(primary)

length(brca_rna_codes) # 1093 barcodes

# Retrieve patients with available mutational data
brca_mut_codes <- brca_samples %>% 
  filter(assay == "BRCA_Mutation-20160128", sample_types == "01") %>%
  pull(primary)

length(brca_mut_codes) # 988 barcodes

# Retrieve samples with both transcriptomic and mutational information
brca_codes <- intersect(brca_rna_codes, brca_mut_codes)
length(brca_codes) # 974 patients

# ==========================================================

#  __  __       _        _   _                   _       _       _        
# |  \/  |_   _| |_ __ _| |_(_) ___  _ __   __ _| |   __| | __ _| |_ __ _ 
# | |\/| | | | | __/ _` | __| |/ _ \| '_ \ / _` | |  / _` |/ _` | __/ _` |
# | |  | | |_| | || (_| | |_| | (_) | | | | (_| | | | (_| | (_| | || (_| |
# |_|  |_|\__,_|\__\__,_|\__|_|\___/|_| |_|\__,_|_|  \__,_|\__,_|\__\__,_|

# Divide patients whether they have a mutation in BRCA1 or BRCA2 (BRCA-positive)

# Retrieve mutational data
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  barcode = brca_codes # Patients with both transcriptomic and mutational data
)

GDCdownload(query) # Run this line once, make sure you are working on the correct directory

maf <- GDCprepare(query)

maf # A tibble with mutational information reported in patients

# Retrieve barcodes from patients with mutations in BRCA1 or BRCA2
brcapos_codes <- unique(maf %>% filter(Hugo_Symbol == "BRCA1" | Hugo_Symbol == "BRCA2") %>% 
                      pull(Tumor_Sample_Barcode))

length(brcapos_codes) # 42 patients

# In this case patient barcodes correpond to mutational data, not clinical data
# To retrieve clinical data we must shorten mutational barcodes to general barcodes 
# Clinical data is attached to the general barcodes
# 
# To read more about TCGA barcodes, enter the following website:
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

brcapos_codes <- substring(brcapos_codes, 1, 12)

# Retrieve patients BRCA-negative (No mutations in BRCA1 or BRCA2)
brcaneg_codes <- setdiff(brca_codes, brcapos_codes)

length(brcaneg_codes) # 932 patients

# ==========================================================

#   ____ _ _       _           _       _       _        
#  / ___| (_)_ __ (_) ___ __ _| |   __| | __ _| |_ __ _ 
# | |   | | | '_ \| |/ __/ _` | |  / _` |/ _` | __/ _` |
# | |___| | | | | | | (_| (_| | | | (_| | (_| | || (_| |
#  \____|_|_|_| |_|_|\___\__,_|_|  \__,_|\__,_|\__\__,_|

metadata_tcga <- as.data.frame(colData(brca)) # Convert clinical data as a data frame
# write.csv(metadata_tcga, "nonfiltered_metadata.csv", row.names = TRUE)
# saveRDS(metadata_tcga, "nonfiltered_metadata.RDS")

# BRCA-positive; pre-menopause

brcapos_menopre <- metadata_tcga %>% filter(patientID %in% brcapos_codes,
                         patient.menopause_status == "pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement)")

nrow(brcapos_menopre) # 11 patients

# BRCA-positive; post-menopause

brcapos_menopos <- metadata_tcga %>% filter(patientID %in% brcapos_codes,
                                            patient.menopause_status == "post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy)")

nrow(brcapos_menopos) # 24 patients

# BRCA-negative; pre-menopause

brcaneg_menopre <- metadata_tcga %>% filter(patientID %in% brcaneg_codes,
                                            patient.menopause_status == "pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement)")

nrow(brcaneg_menopre) # 197 patients

# BRCA-negative; post-menopause

brcaneg_menopos <- metadata_tcga %>% filter(patientID %in% brcaneg_codes,
                                            patient.menopause_status == "post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy)")

nrow(brcaneg_menopos) # 610 patients 

# ==
# Generate a clinical data frame
# ==

# Order barcodes
primary_codes <- c(brcapos_menopre$patientID,
                   brcapos_menopos$patientID,
                   brcaneg_menopre$patientID,
                   brcaneg_menopos$patientID)

# Mutation barcodes
mut_barcodes <- brca_samples %>% filter(assay == "BRCA_Mutation-20160128", sample_types == "01") %>% 
  distinct(primary, .keep_all = TRUE) %>% filter(primary %in% primary_codes) %>% 
  arrange(match(primary, primary_codes)) %>% pull(colname)

# Expression barcodes
expr_barcodes <- brca_samples %>% filter(assay == "BRCA_RNASeq2GeneNorm-20160128", sample_types == "01") %>%
  distinct(primary, .keep_all = TRUE) %>% filter(primary %in% primary_codes) %>%
  arrange(match(primary, primary_codes)) %>% pull(colname)

#BRCA mutation (BRCA-pos patients)
brcapos_mutations_df <- maf %>% filter(Hugo_Symbol == "BRCA1" | Hugo_Symbol == "BRCA2")
brcapos_mutations_df$primary_codes <- substring(brcapos_mutations$Tumor_Sample_Barcode, 1, 12)

brcapos_mutations <- c() 
# Empty vector for BRCA mutations with the order of the BRCA-positive patients
# in "primary_codes" vector

for (barcode in primary_codes[1:35]) { # 35 BRCA-pos patients
  barcode_mut <- brcapos_mutations_df %>% filter(primary_codes == barcode)
  
  mutations <- unique(barcode_mut$Hugo_Symbol)
  
  if ("BRCA2" %in% mutations & "BRCA1" %in% mutations) {
    mutations <- "BRCA1 / BRCA2"
  }
  
  brcapos_mutations <- c(brcapos_mutations, mutations)
}

brcapos_mutations # 35

# Menopause status
menopause <- metadata_tcga %>% filter(patientID %in% primary_codes) %>%
  pull(patient.menopause_status)

menopause <- substring(menopause, 1, 3) # Only "pre" or "pos"

# Age status
age <- metadata_tcga %>% filter(patientID %in% primary_codes) %>%
  pull(Age.at.Initial.Pathologic.Diagnosis)

# Create clinical data frame
clinical_metadata <- data.frame(patienpatient_idclinical_metadata <- data.frame(patients = primary_codes,
                                mutational_barcodes = mut_barcodes,
                                expression_barcodes = expr_barcodes,
                                age_at_initial_diagnosis = age,
                                menopause_status = c(menopause),
                                brca_mutations = c(brcapos_mutations, rep(NA, times = 807)),
                                brca_status = c(rep("+", times = 35), rep("-", times = 807)),
                                experimental_design = c(rep("BRCA-positive; Pre-Menopause", times = 11),
                                                        rep("BRCA-positive; Post-Menopause", times = 24),
                                                        rep("BRCA-negative; Pre-Menopause", times = 197),
                                                        rep("BRCA-positive; Post-Menopause", times = 610))))

write.csv(clinical_metadata, "clinical_metadata.csv", row.names = TRUE)
saveRDS(metadata_tcga, "clinical_metadata.RDS")

# ==========================================================

#  _____                              _                   _       _        
# | ____|_  ___ __  _ __ ___  ___ ___(_) ___  _ __     __| | __ _| |_ __ _ 
# |  _| \ \/ / '_ \| '__/ _ \/ __/ __| |/ _ \| '_ \   / _` |/ _` | __/ _` |
# | |___ >  <| |_) | | |  __/\__ \__ \ | (_) | | | | | (_| | (_| | || (_| |
# |_____/_/\_\ .__/|_|  \___||___/___/_|\___/|_| |_|  \__,_|\__,_|\__\__,_|
#            |_|                                                           

# Retrieve RNA-Seq data from expression barcodes
counts.query <- GDCquery(project = "TCGA-BRCA",
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification", 
                         experimental.strategy = "RNA-Seq",
                         workflow.type = "STAR - Counts",
                         barcode = clinical_metadata$expression_barcodes) # Patients from our experimental design

GDCdownload(counts.query) # Download once, make sure to work on the correct directory

raw.counts_brca <- GDCprepare(query = counts.query)

counts_brca <- raw.counts_brca@assays@data@listData$unstranded # Raw matrix

# Label columns using expression barcodes
colnames(counts_brca) <- colnames(raw.counts_brca)

# Label rows using ENSEMBL
ensg <- unlist(strsplit(rownames(raw.counts_brca), split = "[.]")) # Remove human genome version
ensg <- ensg[c(TRUE, FALSE)]
rownames(counts_brca) <- ensg 

dim(counts_brca) # [1] 60660   842

write.csv(counts_brca, "expression_matrix.csv", row.names = TRUE)
saveRDS(metadata_tcga, "expression_matrix.RDS")

# save.image("~/Documents/BrCanProject/metadata_image.RData")
