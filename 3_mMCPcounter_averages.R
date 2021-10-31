
######################################################################
# __________ mMCP-COUNTER __________ #
# https://github.com/cit-bioinfo/mMCP-counter
######################################################################

#rm(list = ls()) #clear environment
getwd() # get working directory
dev.off() # To reset heatmap render error

library(tidyverse)
library(mMCPcounter)
library(ComplexHeatmap)
library(RColorBrewer)
dir.create('mMCPcounter')

### *** on next use need to insert toggle for average calculation *** ###


######################################################################
# 1) IMPORT -> mMCP OUTPUT
######################################################################

# Import from DESeq2 normalised counts, needs to be matrix and +1 log2 transformed
# Run mMCP-count.estimate function where Gene.Symbol = MGI symbols present as matrix rownames
mMCP.input <- as.matrix(read.csv("DESeq_export/normalised_counts.csv", header=T, row.names=1, sep=","))
mMCP.input <- log2(mMCP.input+1)
mMCP.output <- mMCPcounter.estimate(mMCP.input, features = c("Gene.Symbol")[1])
head(mMCP.output)


######################################################################
# 2) Z-SCORE ACROSS ALL ROWS: CREATE EMPTY MATRIX FOR AVERAGES & POPULATE
######################################################################

# Use gsub to convert column names to unique names shared by treatment {treatment_timepoint}
# Create empty matrix using length of rownames and unique column names in column.list
# Populate empty matrix rownames and column names
print(mMCP.column.list <- unique(gsub(pattern = 'EPA01_[A-Z][0-9]_MOC2_', replacement = '', colnames(mMCP.output))))
mMCP.averages <- matrix(nrow = length(rownames(mMCP.output)), ncol = length(mMCP.column.list))
rownames(mMCP.averages) <- rownames(mMCP.output)
colnames(mMCP.averages) <- mMCP.column.list
mMCP.averages # empty matrix

# Filters by {treatment_timepoint} name,  calculates rowmeans and populates empty matrix with averages
# NOTE: Cannot be gaps in data or the loop will exit with an error, must fill with NA instead
for (a in mMCP.column.list) {
  mMCP.averages[,paste0(a)] <- rowMeans(mMCP.output[ , grepl(paste0(a), colnames(mMCP.output))], na.rm = TRUE)
}
rm('a', 'mMCP.column.list') # remove loop related lists


######################################################################
# 3) Z-SCORE ACROSS ALL ROWS: HEATMAP
######################################################################

# [MANUAL] Reorder the columns and rows if necessary
mMCP.averages.tidy <- mMCP.averages
colnames(mMCP.averages.tidy)
mMCP.averages.tidy <- mMCP.averages.tidy[ ,c('Vehicle_D3',
                                             'ATRi_D3',
                                             'RT_D3',
                                             'ATRiRT_D3',
                                             'Vehicle_D10',
                                             'ATRi_D10',
                                             'RT_D10',
                                             'ATRiRT_D10'
                                             )]

rownames(mMCP.averages.tidy)
mMCP.averages.tidy <- mMCP.averages.tidy[c(3,2,1,4:15),] # moving NK cells (3) to top and CD8 to second (2)
mMCP.averages.tidy <- mMCP.averages.tidy[c(1:11),] # trim vessels, lymphatics, endothelial cells, fibro
rownames(mMCP.averages.tidy)[rownames(mMCP.averages.tidy) == 'Monocytes / macrophages'] <- 'Monocytes/Macrop.'
mMCP.averages.tidy

# [MANUAL] Rename the columns if necessary
colnames(mMCP.averages.tidy) <- c('Vehicle',
                                  'ATRi',
                                  'RT',
                                  'ATRi/RT',
                                  'Vehicle',
                                  'ATRi',
                                  'RT',
                                  'ATRi/RT'
                                  )

# Timepoint Annotation
print(mMCP.annotation <- data.frame(Timepoint =c('Day 3','Day 3','Day 3','Day 3','Day 10','Day 10','Day 10','Day 10')))
heatmap.top.annotation = HeatmapAnnotation('Timepoint' = mMCP.annotation$Timepoint,
                                           col = list('Timepoint' = c('Day 3' = '#fdb361', 'Day 10' = '#00A087FF')),
                                           annotation_name_side = "left",
                                           annotation_name_gp = gpar(fontsize = 8),
                                           simple_anno_size = unit(0.3, "cm"),
                                           gap = unit(0, "cm"),
                                           annotation_legend_param = list('Timepoint' = list( # the legend can be any group
                                                title = "Timepoint",
                                                title_gp = gpar(fontsize = 8, font = 2),
                                                labels_gp = gpar(fontsize = 8),
                                                at = c("Day 3", "Day 10"), 
                                                labels = c("Day 3", "Day 10") # To change the names of categories
                                              )))


# scale
# heatmap
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
                     col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                     heatmap_legend_param = list(
                       title = 'Z-Score',
                       title_gp = gpar(fontsize = 8, font = 2),
                       labels_gp = gpar(fontsize = 8),
                       grid_width = unit(.4, "cm")
                     ),
                     column_split = factor(rep(c("Day3", "Day10"), each = 4), levels = c("Day3", "Day10")),
                     column_title = NULL,
                     border = TRUE,
                     top_annotation = heatmap.top.annotation
)

draw(heatmap.plot, merge_legend = TRUE)


pdf(paste0('mMCPcounter/', 'newer_script', '.pdf'), width=3, height=2.1)
draw(heatmap.plot, merge_legend = TRUE)
dev.off()



# Tidy up at end
#rm(list = ls(pattern = 'mMCP|heatmap'))

