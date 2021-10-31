###################################################
# ----- DIFFERENTIAL EXPRESSION WITH DESeq2 ----- #
###################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # sets working directory based on script location - should be in analysis folder
getwd() # get working directory

library(tidyverse)
library(DESeq2)
library(crayon)

# Note: this is not set up to perform batch correction

#########################################
# STEP 1: DETAILS ON INPUTS AND OUTPUTS #
#########################################

##### INPUT FILES (IN PROJECT ROOT DIRECTORY) ##### 
# countData = "counts_gene.csv" output from prepDE.py3
# colData = sample_names.csv
    # sample_names            = these are the new samples names provided during the alignment script 1 renaming
    # old_sample_name         = if a previous sample name used 
    # treatment_wo_timepoint  = vehicle/PD1/RT/RT_PD1/etc
    # timepoint               = 14d/21d, injected/contralateral, parental/KO, etc
    # treatment               = 14d_vehicle/14d_PD1/21d_vehicle/21d_PD1/etc
# Note: colData$sample_names must be in exactly the same order as the column names containing samples names in countData (counts_gene.csv)


##### INPUT VARIABLES ##### 
control_name='Vehicle' # This is case sensitive (Vehicle, vehicle, Control, control, etc)
timepoints=c('D3', 'D10') # 'timepoints' = 1) timepoint 2) biflank 3) KO/UV/OE/etc 4) different models if single timepoint each
count_filter_number = 24 # Set to filter rows that have less than 1 read per sample - ie, how many samples do you have?

##### OUTPUTS #####
# dds data object for each timepoint
# DESeq_export/normalised_counts.csv
# DESeq_export/timepoint1/ all DGE comparison
# DESeq_export/timepoint2/ all DGE comparison
# DESeq_export/DESeq_output.RData save of environment


########################
# STEP 2: DATA IMPORT #
########################

colData <- read.csv("./../sample_names.csv", row.names="sample_names")
countData <- read.csv("./../counts_gene.csv")

# DESeq2 requires colData columns to be factors
colData$treatment <- factor(colData$treatment, levels = unique(colData$treatment))
colData$timepoint <- factor(colData$timepoint, levels = unique(colData$timepoint))
colData$treatment_wo_timepoint <- factor(colData$treatment_wo_timepoint, levels = unique(colData$treatment_wo_timepoint))

# Stringtie/prepDE.py3 exports the gene_id in the format "ENSEMBL-ID|mgi_symbol"
# This leaves only the mgi_symbol, checks for duplicates, sums duplicates, performs a check, then converts to a matrix
countData$gene_id <- gsub("[A-Z0-9]*\\|", "", countData$gene_id) # find and replace with grep-substitute
countData <- countData[order(countData$gene_id),] # order for duplicate check
countData_qc_duplicates <- countData %>% group_by(gene_id) %>% filter(n() > 1) # qc export of duplicates found
countData <- aggregate(. ~ gene_id, data=countData, FUN=sum) # sums column values when merging duplicates
countData %>% group_by(gene_id) %>% filter(n() > 1) # duplicate check after sum, should be just headers
countData <- countData %>% column_to_rownames(var = "gene_id") %>% as.matrix() # add rownames, convert to matrix

# Additional checks to determine sample names match between columns (countData) and rows (colData)
print(all(rownames(colData) %in% colnames(countData))) # do all row names in sample data (colData) appear in column names for counts (countData)?
print(all(rownames(colData) == colnames(countData))) # are they in the same order? # They need to be
# countData <- countData[, rownames(colData)] # If needed, reorders countData columns to match sample order in colData file - run two lines above again


##################################################
# STEP 3: DESEQ2 FOR LOOP BY TIMEPOINT/FLANK/ETC #
##################################################

for (timepoint in timepoints) {
  
# [1] Create a DESeqDataSet (dds) from count matrix (countData) and sample info (colData) using 
# Note: DESeqDataSetFromMatrix() is specific for ballgown>prepDE.py3 output, needs changed if alignment switched to HTSeq/Salmon/other
# Removes rows with average of 1 read per sample
# Manual: "more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function."
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ treatment)
cat(magenta(paste('Rows Pre-Filtering:',nrow(dds),'\n'))) # genes pre filtering
keep <- rowSums(counts(dds)) >= count_filter_number # results in a TRUE/FALSE output # count_filter_number is at very start and = 1 read per sample
dds <- dds[keep,] # only keeps rows = TRUE
cat(magenta(paste('Rows Post-Filtering:', nrow(dds),'\n'))) # genes post filtering
rm(keep)

# [2] Set "Factor" Levels & DGE Analysis
# 'reference_level' sets reference level for DGE comparison, changes based on timepoint in for loop using paste to change the timepoint on each loop
# Run DGE analysis using DESeq() set to use multiple cores
reference_level=paste0(control_name, '_', timepoint)
dds$treatment <- relevel(dds$treatment, ref = reference_level)
dds <- DESeq(dds, parallel = TRUE)
cat(magenta(paste('List of comparisons at:', timepoint, '\n')))

# [3] Export
dir.create("DESeq_export/", showWarnings = FALSE)
# This code runs every loop, no issues with overwrite as normalised counts always the same irrespective of progression through loops
write.csv(as.data.frame(counts(dds, normalized=T)), file="DESeq_export/normalised_counts.csv")
# all DGE comparisons using for loop to export csv files
directory=paste0('DESeq_export/', timepoint, '/')
dir.create(directory, showWarnings = FALSE)
DGE_comparisons=resultsNames(dds) # lookup results comparisons list for export loop
print(DGE_comparisons)
for (i in DGE_comparisons) {write.csv(as.data.frame(results(dds, name=i)), paste0(directory, i, ".csv"))}

# [4] Rename dds by timepoint
# this renames the dds file based on the timepoint so that it can be used to import results in subsequent analyses
dds_name=paste0('dds_', timepoint)
assign(dds_name, dds)
rm(dds)
cat(green(paste0('DEseq2: ', timepoint, ' Loop Completed\n')))

} # end of timepoints for loop


#####################################################
# STEP 4: Remove Unneeded Data and Save .RDATA file #
#####################################################

rm(dds_name, DGE_comparisons, directory, i, reference_level, countData, countData_qc_duplicates, count_filter_number)
save.image(file = "DESeq_export/DESeq_output.RData")

# End of Script
