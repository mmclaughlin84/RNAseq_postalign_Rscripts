######################################################
# ____________ DESeq2 ___________#
#####################################################

# Manu Colour Scheme
# Control = #808080
# ATRi = #b855d8
# RT = #08bf01
# ATRi/RT = #5657f9

# grey = #e9eef0
# red = #ff4e27
# cyan = #4DBBD5FF
# green = #00A087FF
# light orange = #fdb361

rm(list = ls()) #clear environment

getwd() # get working directory
setwd("/Users/mclaughinm/Desktop/RNAseq_EPMOC2/analysis") # set working directory as *ANALYSIS* folder
list.files()

# All updated for Rv4.03
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)
library(dendextend)
library(crayon)
dev.off() # To reset heatmap render error
# Future use: # batch correction with ComBat #normalisation with fRMA

###########################################
### STEP 1: DETAILS ON DATA PREP ###
###########################################

##### INPUTS (IN PROJECT ROOT DIRECTORY) ##### 

# countData = "counts_ENonly_gene.csv" output from prepDE
# colData = sample_names.csv
    # sample_names           = {experimnetid}_{model}_{timepoint}_{treatment}_{replicate} # ONLY _ do not use hyphens
    # old_sample_name         = if a previous sample name used 
    # treatment_wo_timepoint  = vehicle/PD1/RT/RT_PD1/etc
    # timepoint               = 14d/21d, injected/contralateral, parental/KO, etc
    # treatment               = 14d_vehicle/14d_PD1/21d_vehicle/21d_PD1/etc

##### INPUT CONTROL & TIMEPOINT VARIABLES ##### 

control_name='Vehicle'
timepoints=c('D3', 'D10') # 'timepoints' = 1) timepoint 2) biflank 3) KO/UV/OE/etc
count_filter_number = 24 # filter rows that have less than 1 read per sample

##### OUTPUT #####
# dds data object for each timepoint
# DESeq_export/normalised_counts.csv
# DESeq_export_timepoint1/ all DGE comparison
# DESeq_export_timepoint2/ all DGE comparison
##### END ##### 




###########################################
### STEP 2: DATA IMPORT ###
###########################################

### 1) Import
# DESeqDataSetFromMatrix() for ballgown>prepDE.py, different for HTSeq/Salmon/Sailfish/kallisto/RSEM
# countData = tidied count data with gene_ids in column 1 used as row.names for matrix
# sampleData = 1) sample_names as 1st column must matching first row in countData
# sampleData = 1) treatment is groups to be used for DGE analysis
# Column headers cannot start with a number, or contain - 
colData <- read.csv("./../sample_names.csv", row.names="sample_names")
countData <- read.csv("./../counts_ENonly_gene.csv")

### 1b) Name Crosschecking - additional checks require conversion to matrix first
colnames(countData) <- gsub("\\.", "_", colnames(countData)) # Convert . and - to _ it makes life easier # REMEMBER, MUST MATCH colData!!!!!!!!!
rownames(colData) <- gsub("-", "_", rownames(colData)) # Make sample names in colData match countData

### Convert colData columsn to factors
str(colData)
colData$treatment <- factor(colData$treatment, levels = unique(colData$treatment))
colData$timepoint <- factor(colData$timepoint, levels = unique(colData$timepoint))
colData$treatment_wo_timepoint <- factor(colData$treatment_wo_timepoint, levels = unique(colData$treatment_wo_timepoint))

# Remove ENSEMBL-ID| and duplicates, convert to matrix
countData$gene_id <- gsub(".*\\|", "", countData$gene_id) # find and replace with grep-substitute
countData <- countData[order(countData$gene_id),] # order for duplicate check
countData_qc_duplicates <- countData %>% group_by(gene_id) %>% filter(n() > 1) # qc export of duplicates found
countData <- aggregate(. ~ gene_id, data=countData, FUN=sum) # sums column values when merging duplicates
countData %>% group_by(gene_id) %>% filter(n() > 1) # duplicate check after sum, should be just headers
countData <- countData %>% column_to_rownames(var = "gene_id") %>% as.matrix() # add rownames, convert to matrix

# Additional checks to determine sample names match between columns (countData) and rows (colData)
print(all(rownames(colData) %in% colnames(countData))) # do all row names in sample data (colData) appear in column names for counts (countData)?
print(all(rownames(colData) == colnames(countData))) # are they in the same order? # They need to be
# countData <- countData[, rownames(colData)] # This will reorder the countData columns to match the sample order provided in the colData file
# all(rownames(colData) == colnames(countData)) # It reorders sample names and numbers in the column, this has been triple checked, but do it every time to be certain





######################################################
### STEP 3: DESEQ2 FOR LOOP BY TIMEPOINT/FLANK/ETC ###
######################################################

for (timepoint in timepoints) {
  
### 2) DEseq matrix and filtering
# Create a DESeqDataSet from count matrix (countData) and sample labels (colData)
# Removes rows with average of 1 read per sample ("more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ treatment)
cat(magenta(paste('Rows Pre-Filtering:',nrow(dds),'\n'))) # no. pre filtering
keep <- rowSums(counts(dds)) >= count_filter_number # results in a TRUE/FALSE output # count_filter_number is at very start and = 1 read per sample
dds <- dds[keep,] # only keeps rows = TRUE
cat(magenta(paste('Rows Post-Filtering:', nrow(dds),'\n'))) # no. post filtering
rm(keep)


### 3) Set "Factor" Levels & DGE Analysis
# sets reference level for DGE comparison, changes based on timepoint in for loop
# way to add second variable, such as batch (for batch correction?)
# Run DGE analysis with multiple cores
reference_level=paste0(control_name, '_', timepoint)
dds$treatment <- relevel(dds$treatment, ref = reference_level)
dds <- DESeq(dds, parallel = TRUE)
cat(magenta(paste('List of comparisons at:', timepoint, '\n')))



### 4) Export
# a) Normalised counts 
# these two lines will duplicate on each loop, but no issues with overwrite as normalised counts always the same irrespective of timepoint and treatment setup factor/levels setup
dir.create("DESeq_export/", showWarnings = FALSE) 
write.csv(as.data.frame(counts(dds, normalized=T)), file="DESeq_export/normalised_counts.csv") # normalised not be impacted by timepoint setting
#b) all DGE comparisons using for loop to export csv files
directory=paste0('DESeq_export/', timepoint, '/')
dir.create(directory, showWarnings = FALSE)
DGE_comparisons=resultsNames(dds) # lookup results comparisons for export loop
print(DGE_comparisons)
for (i in DGE_comparisons) {
  write.csv(as.data.frame(results(dds, name=i)), paste0(directory, i, ".csv"))
}

### 5) Rename dds by timepoint
# this renames the dds file based on the timepoint so that it can be used to import results in subsequent analyses
print('Note: Normalised counts are the same irrespective of dds timepoint used')
dds_name=paste0('dds_', timepoint)
assign(dds_name, dds)
rm(dds)

cat(green(paste0('DEseq2: ', timepoint, ' Loop Completed\n')))

}



###############################################################
### STEP 4: Removal of Unneeded Data and Saving .RDATA file ###
###############################################################

rm(dds_name, DGE_comparisons, directory, i, reference_level, countData, countData_qc_duplicates, count_filter_number)
save.image(file = "DESeq_export/DESeq_output.RData")
#load(file = "DESeq_export/DESeq_output.RData")

###############################
# _______At This Point_______ #
###############################
### >>>>>> PCA script + All heatmap + DEG count <<<<< ###
### >>>>> GSEA/TopGO/GO# scripts <<<<< ###
### >>>>> Volcano Plot <<<<< ### not yet written, copy VR
### >>>>> Heatmaps: GO term based <<<<< ####


