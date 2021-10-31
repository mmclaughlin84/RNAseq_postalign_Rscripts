####################################################################
#_____HEATMAPS: log2 z-scaled values based on custom gene list_____#
####################################################################

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize) #colorRamp2
library(DESeq2)
suppressWarnings(dir.create('Heatmaps'))
setwd("/Users/mclaughinm/Desktop/RNAseq_MP106_RT_STINGa/analysis") # set working directory as *ANALYSIS* folder
rm(list = ls()) #clear environment
load(file = "DESeq_export/DESeq_output.RData") # Avoids rerunning DESeq2


######################################################################
# 0) Variables
######################################################################

var_panel_name='general'
var_what_timepoint="21d" # Select timepoint/biflank id(s) in nomenclature
var_averages_or_not='all_values' #averages|all_values, used for if function and export file name
var_regex_start='^everything_before_treatment_time' # Regular expression to leave only {treatment}_{timepoint}
var_regex_end='_r[0-9]$|_rep[0-9]$' # trim old _r[0-9] naming system if needed
# [1st list] = order samples are actually in. [2nd list (levels)] = order they SHOULD be in
var_simple_treatment_names=factor(c('Vehicle', 'BI880', 'RT', 'RT-BI880'), levels = c('Vehicle', 'BI880', 'RT', 'RT-BI880'))


######################################################################
# 1) Import Custom Gene List
######################################################################

# Load custom gene list & create .long version as required for certain steps
# Print duplicate entry check (will be blank if no duplicate gene ids)
c.gene.list <- read_csv(file = paste0('Heatmaps/custom_gene_list_', var_panel_name, '.csv'))
c.gene.list.long <- c.gene.list %>% t() %>% data.frame()
c.gene.list.long$category <- rownames(c.gene.list.long)
c.gene.list.long <- pivot_longer(c.gene.list.long, cols = 1:(ncol(c.gene.list.long)-1), values_to = c('gene_id'), values_drop_na = TRUE, names_to = NULL)
print(c.gene.list.long[ duplicated(c.gene.list.long$gene_id), ])


######################################################################
# 2) Load Scaled Normalised Counts
######################################################################

# Import normalised counts and filter by timepoint
# Log2+1 values and convert to matrix
c.counts.norm <- data.frame(read.csv('DESeq_export/normalised_counts.csv', row.names = 1)) %>% tibble::rownames_to_column("gene_id")
c.counts.norm <- dplyr::select(c.counts.norm, matches(paste0(var_what_timepoint,'|gene_id'))) # Select for different timepoints or multiple timepoints
rownames(c.counts.norm) <- c.counts.norm$gene_id # add rownames
c.counts.norm <- c.counts.norm[ ,-1] # delete gene_id column, [subset] does not allow -'gene_id'
c.counts.norm <- log2(as.matrix(c.counts.norm)+1) # log2+1 and convert to matrix
colnames(c.counts.norm) <-gsub(pattern = var_regex_end, replacement = '', colnames(c.counts.norm)) # remove _r/rep# at end
colnames(c.counts.norm) <-gsub(pattern = var_regex_start, replacement = '', colnames(c.counts.norm)) # leave only {treatment}_{timepoint}


######################################################################
# 3) CALCULATE AVERAGES OF NORMALISED COUNTS IF VARIABLE SELECTED
######################################################################

# If var_averages_or_not == averages then c.counts.norm will be overwritten by c.counts.norm.averages intermediate
# gsub has earlier converted column names to {treatment_timepoint} only based on var_regexs
# Create empty matrix using length of rownames and unique column names in column.list and copying row/column names
# For loop filters by {treatment_timepoint} name,  calculates rowmeans and populates empty matrix with averages
# NOTE: Cannot be gaps in data or the loop will exit with an error, must fill with NA instead

if (var_averages_or_not == 'averages') {
  print('Unique Sample Names For Average Calculation Are:')
  print(column.list <- unique(colnames(c.counts.norm)))
  c.counts.norm.averages <- matrix(nrow = length(rownames(c.counts.norm)), ncol = length(column.list))
  rownames(c.counts.norm.averages) <- rownames(c.counts.norm)
  colnames(c.counts.norm.averages) <- column.list
  for (a in column.list) {
    c.counts.norm.averages[,paste0(a)] <- rowMeans(c.counts.norm[ , grepl(paste0(a), colnames(c.counts.norm))], na.rm = TRUE)
  }
  c.counts.norm <- c.counts.norm.averages
  print('NOTE: Calculation of Averages Requested and Completed')
} else {print('WARNING: Calculation of Averages of Normalised Counts Not Requested')}
suppressWarnings(rm('a', 'column.list', 'c.counts.norm.averages', 'var_regex_start', 'var_regex_end')) # tidy up


######################################################################
# 4) Load DEG Associated Adjusted P-Values
######################################################################

### Assumes 4 groups: could convert to for loop but leaving like this for flexibility ###
# From var_what_timepoint extracts correct comparison for *only* that timepoint versus matching vehicle
# Merged adjusted p-values into c.adjpvalues and renames to var_simple_treatment_names minus 1st vehicle/control name
c.DEG.list <- resultsNames(get(paste0('dds_',var_what_timepoint)))[grep(paste0(var_what_timepoint,'_vs_'), resultsNames(get(paste0('dds_',var_what_timepoint))))] # Pulls timepoint
c.temp1 <- as.data.frame(results(get(paste0('dds_',var_what_timepoint)), name=c.DEG.list[1])) %>% tibble::rownames_to_column("gene_id") # best df <- matrix
c.temp2 <- as.data.frame(results(get(paste0('dds_',var_what_timepoint)), name=c.DEG.list[2])) %>% tibble::rownames_to_column("gene_id") # best df <- matrix
c.temp3 <- as.data.frame(results(get(paste0('dds_',var_what_timepoint)), name=c.DEG.list[3])) %>% tibble::rownames_to_column("gene_id") # best df <- matrix
suppressWarnings(rm(c.adjpvalues)) # thought there might have been a clash when rerunning code so removing if exists just incase
c.adjpvalues <-  c.temp1["gene_id"]
c.adjpvalues[c.DEG.list[1]] <- c.temp1$padj
c.adjpvalues[c.DEG.list[2]] <- c.temp2$padj
c.adjpvalues[c.DEG.list[3]] <- c.temp3$padj
rm(list = ls(pattern = 'c.temp'), 'c.DEG.list')

  
######################################################################
# 4) Filtering Scaled Gene Expression using Custom Gene List
######################################################################

# Filter normalised counts so only those gene_ids in the external custom gene list provided remain
# What genes from the RNAseq are NOT in the custom gene list?
# Filter the custom gene to show those NOT in the RNAseq data 
c.FILT.norm <- c.counts.norm[rownames(c.counts.norm) %in% c.gene.list.long$gene_id, ]
not.in.RNAseq <- !(c.gene.list.long$gene_id %in% rownames(c.FILT.norm)) 
print(not.in.RNAseq <- c.gene.list.long$gene_id[not.in.RNAseq] ) # NOT present in the RNAseq dataset: NEED REMOVED

if (length(not.in.RNAseq) == 0) { print('All custom genes present in RNAseq data')
  } else {
# gene_ids not in RNAseq are collapsed into a long (or) string and lapplied across dataframe using gsub
# rows with all NAs removed,the original custom gene list is 1) overwritten, 2) made long, 3) filter check repeated
not.in.RNAseq.collapsed <- paste(unlist(not.in.RNAseq), collapse = "|")
c.gene.list <- as.data.frame(lapply(c.gene.list, function(x) {gsub(pattern = not.in.RNAseq.collapsed, replacement = NA, x)}))
c.gene.list <- c.gene.list[rowSums(is.na(c.gene.list)) != ncol(c.gene.list), ]
rm('not.in.RNAseq','not.in.RNAseq.collapsed')

# repeating pivot wide step without duplicate check
c.gene.list.long <- c.gene.list %>% t() %>% data.frame()
c.gene.list.long$category <- rownames(c.gene.list.long)
c.gene.list.long <- pivot_longer(c.gene.list.long, cols = 1:(ncol(c.gene.list.long)-1), values_to = c('gene_id'), values_drop_na = TRUE, names_to = NULL)
}

# Repeating check step, anything left is in custom gene list but ***NOT MGI*** list and needs ***PERMENANTLY REMOVED***
c.FILT.norm <- c.counts.norm[rownames(c.counts.norm) %in% c.gene.list.long$gene_id, ]
not.in.RNAseq <- !(c.gene.list.long$gene_id %in% rownames(c.FILT.norm)) 
print(not.in.RNAseq <- c.gene.list.long$gene_id[not.in.RNAseq] )

# These two lists need to be the same
print(rownames(c.FILT.norm)) 
print(sort(c.gene.list.long$gene_id)) 

# Having carried out the cross check, back to the this hamfisted reordering step that requires a double t()
### *** these are already in order, coincidence or needed?
c.FILT.norm <- t(c.FILT.norm)
c.FILT.norm <- c.FILT.norm[, c.gene.list.long$gene_id] # This will reorder the countData columns to match the sample order provided in the colData file
c.FILT.norm <- t(c.FILT.norm)

# Sanity checks on order
all(rownames(c.FILT.norm) == c.gene.list.long$gene_id) # FINAL order crosscheck 1) all present 2) and in matching order
print(rownames(c.FILT.norm)) # These two lists need to be the same
print(c.gene.list.long$gene_id) # These two lists need to be the same


######################################################################
# 5) Filter DGE matrix by Custom Gene List
######################################################################

# The main cross between the gene_ids present in the RNAseq dataset and the custom gene list is above
# Not necessary to repeat, but there is a sanity check at the end  for  gene_id of expression matrix vs the DEG matrix
c.FILT.pvalue <- c.adjpvalues[c.adjpvalues$gene_id %in% c.gene.list.long$gene_id, ]
rownames(c.FILT.pvalue) <- c.FILT.pvalue$gene_id
c.FILT.pvalue <- c.FILT.pvalue[,-1]
c.FILT.pvalue <- t(c.FILT.pvalue)
c.FILT.pvalue <- c.FILT.pvalue[, c.gene.list.long$gene_id] # This will reorder the countData columns to match the sample order provided in the colData file
c.FILT.pvalue <- as.data.frame(t(c.FILT.pvalue))
print(all(rownames(c.FILT.norm) == rownames(c.FILT.pvalue))) # Double checking in the same order? They need to be!

######################################################################
# 6) Renaming/Reordering Sample Names [EXTREME CAUTION AT THIS STEP]
######################################################################

# Skip renaming if averages not requested
if (var_averages_or_not == 'averages') {
  # Filtered Counts: Sample name change
  colnames(c.FILT.norm) # CAUTION: Make sure these are in the correct order before overwrite
  var_simple_treatment_names # CAUTION: Make sure these are in the correct order before overwrite
  colnames(c.FILT.norm) <- var_simple_treatment_names # overwrite
  
  # Filtered Counts: reorder sample names by *factor levels* in var_simple_treatment_names (+print to track)
  print(colnames(c.FILT.norm))
  c.FILT.norm <- c.FILT.norm[ ,levels(var_simple_treatment_names)]
  print(colnames(c.FILT.norm))
  
  # Adjusted P-value: Sample name change
  colnames(c.FILT.pvalue) # CAUTION: As above, make sure these are in the correct order before overwrite
  as.character(var_simple_treatment_names[2:4]) # CAUTION: As above, make sure these are in the correct order before overwrite
  colnames(c.FILT.pvalue) <- var_simple_treatment_names[2:length(var_simple_treatment_names)]
  
  # Adjusted P-value: reorder sample names by *factor levels* in var_simple_treatment_names (+print to track)
  print(colnames(c.FILT.pvalue))
  c.FILT.pvalue <- c.FILT.pvalue[ ,setdiff(levels(var_simple_treatment_names), 'Vehicle')] # Need to exclude Vehicle from list of factors, but assumes this is ALWAYS called 'Vehicle'
  print(colnames(c.FILT.pvalue))
  print('COLUMN RENAMING: Averages requested and columns renamed as listed in var_simple_treatment_names')
} else {print('COLUMN RENAMING: Averages NOT requested, column names not altered')}


######################################################################
# 7a) Annotation: CATEGORIES
######################################################################

# Setup Row Split numbers (manual) and factors as levels from custom list
# Automated Colour Assignment
c.numbers <- unname(colSums(!is.na(c.gene.list))) # generates counts from the columns of the custom gene list
c.factors <- levels(factor(colnames(c.gene.list), levels = colnames(c.gene.list))) # turns colnames into factors
c.category.colours <- data.frame(categories = c.factors, colours = brewer.pal(n = 12, name = 'Paired')[(1+2):(length(c.factors)+2)]) # This works up to ten categories
c.category.colours <- tibble::deframe(c.category.colours) # complex heatmap needs list of vectors not dataframe

# Category annotation code
c.annotation.left = rowAnnotation(category = c.gene.list.long$category,
                                  col = list(category = c.category.colours),
                                  show_annotation_name = c(bar = FALSE),
                                  simple_anno_size = unit(0.1, "cm"),
                                  gap = unit(0, "cm"),
                                  annotation_legend_param = list(
                                    title = 'Category',
                                    title_gp = gpar(fontsize = 7, font = 2),
                                    labels_gp = gpar(fontsize = 7),
                                    at = c(c.factors), # c.factors.inorder???
                                    grid_height = unit(.3, "cm"),
                                    grid_width = unit(.3, "cm")))


######################################################################
# 7b) Heatmap: ADJPVALUES
######################################################################

# convert to matrix
# make anything non-significant an NA value to allow to be grey and not to screw up numeric class
c.FILT.pvalue.matrix <- as.matrix(c.FILT.pvalue)
c.FILT.pvalue.matrix[c.FILT.pvalue.matrix > 0.05] <- NA

draw(c.heatmap.right <- Heatmap(log10(c.FILT.pvalue.matrix),
                           cluster_rows = FALSE, cluster_columns= FALSE, # both turned off
                           column_names_side = "top",
                           row_names_side = 'left',
                           row_names_gp = gpar(fontsize = 7),
                           column_names_gp = gpar(fontsize = 7),
                           col = colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(50), #RdYlBu #Blues(100) # how to clip this for log2
                           width = unit(ncol(c.FILT.pvalue.matrix)*0.3, "cm"),
                           height = unit(nrow(c.FILT.pvalue.matrix)*0.225, "cm"),
                           row_split = factor(rep(c.factors, c.numbers), levels = c.factors),
                           row_title = NULL,
                           row_title_gp = gpar(fill = c(NA), border = NA, fontsize = 7),
                           heatmap_legend_param = list(
                             title = 'log10 adjpvalue',
                             title_gp = gpar(fontsize = 7, font = 2),
                             labels_gp = gpar(fontsize = 7),
                             grid_width = unit(.4, "cm"),
                             legend_height = unit(1.4, "cm")) # This is the minimum height for some reason
))

######################################################################
# 7c) Heatmap: LOG2+1 Z-SCALED COUNTS
######################################################################

# Prints just z-scaled counts followed by draw step with adjusted p-value heatmap added
draw(c.heatmap <- Heatmap(t(scale(t(c.FILT.norm))),
                     cluster_rows = FALSE, # MUST BE TURNED OFF
                     cluster_columns= FALSE, # MUST BE TURNED OFF
                     column_names_side = "top",
                     row_names_side = "left",
                     row_names_gp = gpar(fontsize = 7),
                     column_names_gp = gpar(fontsize = 7),
                     col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50), #RdYlBu #Blues(100) # how to clip this for log2
                     width = unit(ncol(c.FILT.norm)*0.3, "cm"),
                     height = unit(nrow(c.FILT.norm)*0.225, "cm"),
                     row_split = factor(rep(c.factors, c.numbers), levels = c.factors),
                     row_title = NULL,
                     row_title_gp = gpar(fill = c(NA), border = NA, fontsize = 7),
                     heatmap_legend_param = list(
                       title = 'Z-Score',
                       title_gp = gpar(fontsize = 7, font = 2),
                       labels_gp = gpar(fontsize = 7),
                       grid_width = unit(.4, "cm"),
                       legend_height = unit(1.74, "cm")), # This is the minimum height for some reason),
                     left_annotation = c.annotation.left
))

draw(c.heatmap + c.heatmap.right)

######################################################################
# 8) Export
######################################################################

pdf(paste0('Heatmaps/', var_averages_or_not, '_', var_panel_name, '_', var_what_timepoint, '.pdf'), width=4, height=7)
draw(c.heatmap + c.heatmap.right)
dev.off()
