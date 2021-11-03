######################################################################
# __________ mMCP-COUNTER __________ #
# https://github.com/cit-bioinfo/mMCP-counter
######################################################################

# There's a lot of manual code in this related to sample names, sample order, heatmap labeling of timepoints

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # sets working directory based on script location - should be in analysis folder
getwd() # check working 

dev.off() # To reset heatmap render error that sometimes happens

library(tidyverse)
library(mMCPcounter)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize) #colorRamp2
suppressWarnings(dir.create('mMCPcounter'))

######################################################################
# 1) IMPORT -> mMCP OUTPUT
######################################################################

# Import from DESeq2 normalised counts, needs to be matrix and +1 log2 transformed
# Run mMCP-count.estimate function where Gene.Symbol = MGI symbols present as matrix rownames
mMCP.input <- as.matrix(read.csv("DESeq_export/normalised_counts.csv", header=T, row.names=1, sep=","))
mMCP.input <- log2(mMCP.input+1)
mMCP.output <- mMCPcounter.estimate(mMCP.input, features = c("Gene.Symbol")[1])
head(mMCP.output)
write.csv(mMCP.output, file = 'mMCPcounter/mMCP_table_estimates_all.csv')

######################################################################
# 2) Z-SCORE ACROSS ALL ROWS: CREATE EMPTY MATRIX FOR AVERAGES & POPULATE
######################################################################

# Use gsub to convert column names to unique names shared by treatment {treatment_timepoint}
# Create empty matrix using length of rownames and unique column names in column.list
# Populate empty matrix rownames and column names (this method keeps the current column order, unlike using aggregate with fun=mean)
mMCP.column.list <- unique(gsub(pattern = 'JKVR002_4434_[A-Z][0-9]{2}_([a-z]*_[A-Za-z0-9]*)_[0-9]', replacement = '\\1', colnames(mMCP.output)))
mMCP.column.list # replicates no longer have unique names - needed for calculating averages
mMCP.averages <- matrix(nrow = length(rownames(mMCP.output)), ncol = length(mMCP.column.list))
rownames(mMCP.averages) <- rownames(mMCP.output)
colnames(mMCP.averages) <- mMCP.column.list
mMCP.averages # empty matrix

# Filters by {treatment_timepoint} name,  calculates rowmeans and populates empty matrix with averages
# NOTE: Cannot be gaps in data or the loop will exit with an error, must fill with NA instead
for (a in mMCP.column.list) {
  mMCP.averages[,paste0(a)] <- rowMeans(mMCP.output[ , grepl(paste0(a), colnames(mMCP.output))], na.rm = TRUE)
}
mMCP.averages # now full matrix
write.csv(mMCP.averages, file = 'mMCPcounter/mMCP_table_averages.csv')
rm(a, mMCP.column.list, mMCP.input, mMCP.output) # remove loop related lists


######################################################################
# 3) Z-SCORE ACROSS ALL ROWS: HEATMAP
######################################################################

# [MANUAL] Reorder the columns and rows if necessary - but you've followed the instructions so far so shouldn't have to!
mMCP.averages.tidy <- mMCP.averages
colnames(mMCP.averages.tidy)

# Skip this if reorder is not needed - Note you need to reorder before you rename the samples below otherwise they won't be unique
mMCP.averages.tidy <- mMCP.averages.tidy[ , c('injected_control',
                                              'injected_aPD1',
                                              'injected_RP1',
                                              'injected_RP1aPD1',
                                              'contralateral_control',
                                              'contralateral_aPD1',
                                              'contralateral_RP1',
                                              'contralateral_RP1aPD1'
)]

# [MANUAL] Rename the columns removing injected/contralateral as this is inserted in the 
colnames(mMCP.averages.tidy) <- c('Control',
                                  'aPD1',
                                  'RP1',
                                  'RP1-aPD1',
                                  'Control',
                                  'aPD1',
                                  'RP1',
                                  'RP1-aPD1'
                                  )

# [MANUAL] Remove/Select/Reorder which populations to plot
data.frame(population = rownames(mMCP.averages.tidy)) # so you can see what the numbers are
mMCP.averages.tidy <- mMCP.averages.tidy[c(2,1,3,4,5,6,7,8,9), ]
rownames(mMCP.averages.tidy)[rownames(mMCP.averages.tidy) == 'Monocytes / macrophages'] <- 'Monocytes/Macro.' # Shorten a bit
mMCP.averages.tidy
write.csv(mMCP.averages.tidy, file = 'mMCPcounter/mMCP_table_averages_selected.csv')


# [Manual] Timepoint Annotation - quite a few places where the timepoint needs to be changed here
mMCP.annotation <- data.frame(Timepoint =c('Injected','Injected','Injected','Injected',
                                                 'Contralateral','Contralateral','Contralateral','Contralateral'))

heatmap.top.annotation = HeatmapAnnotation('Timepoint' = mMCP.annotation$Timepoint,
                                           col = list('Timepoint' = c('Injected' = '#fdb361', 'Contralateral' = '#00A087FF')),
                                           annotation_name_side = "left",
                                           annotation_name_gp = gpar(fontsize = 8),
                                           simple_anno_size = unit(0.3, "cm"),
                                           gap = unit(0, "cm"),
                                           annotation_legend_param = list('Timepoint' = list( # the legend can be any group
                                                title = "Timepoint",
                                                title_gp = gpar(fontsize = 8, font = 2),
                                                labels_gp = gpar(fontsize = 8),
                                                at = c("Injected", "Contralateral"), 
                                                labels = c("Injected", "Contralateral") # To change the names of categories
                                              )))


# Heatmap Main Body (again timepoint needs changed)
mMCP.averages.scaled <- t(scale(t(mMCP.averages.tidy)))
heatmap.plot <- Heatmap(mMCP.averages.scaled,
                        cluster_rows = FALSE, # MUST BE TURNED OFF
                     cluster_columns = FALSE, # MUST BE TURNED OFF
                     width = unit(12*0.25, "cm"),
                     height = unit(14*0.25, "cm"),
                     column_names_side = "top",
                     row_names_side = "left",
                     column_dend_side = "top",
                     row_dend_side = "right",
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     #col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), # changing to the below colour scheme for cell populations
                     col = colorRamp2(seq(min(mMCP.averages.scaled), max(mMCP.averages.scaled), length = 3), c("blue", "#EEEEEE", "red"), space = "sRGB"),
                     heatmap_legend_param = list(
                       title = 'Z-Score',
                       title_gp = gpar(fontsize = 8, font = 2),
                       labels_gp = gpar(fontsize = 8),
                       grid_width = unit(.4, "cm")
                     ),
                     column_split = factor(rep(c("Injected", "Contralateral"), each = 4), levels = c("Injected", "Contralateral")),
                     column_title = NULL,
                     border = TRUE,
                     top_annotation = heatmap.top.annotation
)

draw(heatmap.plot, merge_legend = TRUE)


pdf(paste0('mMCPcounter/', 'mMCP_counter_plot', '.pdf'), width=3.7, height=2.2)
draw(heatmap.plot, merge_legend = TRUE)
dev.off()

# Clear Environment
rm(list = ls(pattern = 'mMCP|heatmap'))

# end of script