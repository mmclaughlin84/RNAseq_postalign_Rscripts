################################################ 
#_____Basic QC: PCA/Heatmap/#DEGs_____ #
################################################

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # sets working directory based on script location - should be in analysis folder
getwd() # get working directory

# When you created sample names, they were called something like this
# JKVR002_4434_A04_injected_control_4
# This is long, so we want to clip the projectId, model name and sampleID off the start
# Example for above clip='JKVR002_4434_[A-Z][0-9]{2}_'
clip='JKVR002_4434_[A-Z][0-9]{2}_'

##### OUTPUTS #####
# PCA/PCA_plot.pdf      # PCA plot
# PCA/heatmap_all.pdf   # Heatmap all genes
# 


############################
# 1) PCA
############################

# Watch this video on PCA from StatQuest: https://www.youtube.com/watch?v=0Jp4gsfOLMs&t=364s

# Load DESeq2 environment to start
# This imports normalised counts from the csv file
load(file = "DESeq_export/DESeq_output.RData") # Avoids rerunning DESeq2 - is >50Mb so if offsite on VPN can take a moment to load
pcaCounts <- data.frame(read.csv('DESeq_export/normalised_counts.csv', row.names = 1)) %>% tibble::rownames_to_column("gene_id")

# Filter Normalised Counts (rows with zero variance cause an error)
rownames(pcaCounts) <- pcaCounts$gene_id # add rownames from gene id column
pcaCounts <- pcaCounts[ ,-1] # delete gene_id column
pcaCounts <- as.matrix(pcaCounts) # make matrix
keep <- rowSums(pcaCounts) >= (ncol(pcaCounts)*10) # results in a TRUE/FALSE output # set to 10 reads per sample (ie column). Has an impact!
pcaCounts <- pcaCounts[keep,] # only keeps rows = TRUE
rm(keep)

# prcomp (Principal Components Analysis) expects samples as rows and genes as columns so matrix needs to be transposed
# basic quick plot of pc1 and pc2
pca <- prcomp(t(pcaCounts), scale=TRUE)
plot(pca$x[,1], pca$x[,2]) # basic plot

# Extract PrinComp1 and PrinComp2
# add timepoint and treatment from colData to make a properly labeled plot
pca.plot.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1],Y=pca$x[,2]) 
pca.plot.data$timepoint <- colData$timepoint
pca.plot.data$treatment_wo_time <- colData$treatment_wo_time
pca.plot.data # this is what the resulting dataframe looks like

# PCA: calculate how much variation is explained by each principle component (and show on scree plot)
# Shows decreasing variance explained moving through the principle components
pca.var <- pca$sdev^2 # sq of std dev of each pc = variance
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # calc as %
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

# PCA Plot using ggplot
ggplot(data=pca.plot.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(color = treatment_wo_time, shape = timepoint)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")

# Create PCA directory and save PCA plot
dir.create("PCA/", showWarnings = FALSE)
ggsave('PCA/PCA_plot.pdf', plot = last_plot(), width = 15, height = 10, units = 'cm', dpi = 300)

# rename pcaCounts to normalisedCounts for heatmap step below, and cleanup environment using rm()
normalisedCounts <- pcaCounts # 1) matrix 2) filtered - for use in heatmap below
rm(list = ls(pattern = 'pca'))


############################
# 2) Heatmap All Genes
############################

# Uses ComplexHeatmap and RColorBrewer
# Make a copy of normalisedCounts to clip the long experiment details from the sample names using gsub
# transposes, then scales, then transposes back again, as scaling is calculated on columns not rows - and our gene counts are in rows
# scale = - pop. mean / pop. sd, ie zscoring

normalisedCounts_clipped <- normalisedCounts
colnames(normalisedCounts_clipped) <- gsub(clip, '', colnames(normalisedCounts_clipped))

heatmap_plot <- Heatmap(t(scale(t(normalisedCounts_clipped))),
        column_title = 'QC: Heatmap all counts',
        name = "Row Z-Score",
        cluster_rows = TRUE, # TURNED OFF DUE TO LOW NUMBERS
        cluster_columns= TRUE,
        show_column_dend = TRUE,
        row_dend_side = "left",
        column_dend_side = "top",
        column_names_side = "top",
        clustering_distance_rows = "euclidean",
        column_dend_reorder = T,
        row_dend_reorder = T,
        row_names_gp = gpar(fontsize = 1),
        column_names_gp = gpar(fontsize = 7),
        col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100) # bluered(100)
        #heatmap_width = unit(9, "cm"),
        #heatmap_height = unit(6+0.32*nrow(GO_filtered_counts), "cm")
)

# Save heatmap as PDF
pdf('PCA/heatmap_all.pdf', width=5, height=5) # 'open' the pdf saving function
draw(heatmap_plot) # this can take some time as it's drawing a lot of data
dev.off() # tells R you're done drawing and to finish 'saving' the pdf

#cleanup environment
rm(heatmap_plot, normalisedCounts, normalisedCounts_clipped)


############################
# 3) DEGs Count
############################

# 'timepoints' and 'control_name' variables should be loaded into the environment already from .Rdata file at the start
# This is a surprisingly complicated script for what it does and the plots are not all that nice as a 'final' publication graph
# You can take the exported .csv tables and use them to make a custom plot for a figure

for (timepoint in timepoints) {

# Takes dds_timepoint and converts it to just dds to ease looping script
dds_by_timepoint=paste0('dds_', timepoint)
dds='dds'
assign(dds, get(dds_by_timepoint))
rm(dds_by_timepoint)
  
# Get comparison list from resultsNames(dds) and create a blank matrix
(temp.DGE.comparisons=resultsNames(dds)) # lookup results for next line
temp.DGE.count <- matrix(nrow = length(temp.DGE.comparisons), ncol = 2)
rownames(temp.DGE.count) <- temp.DGE.comparisons
colnames(temp.DGE.count) <- c('increase', 'decrease')

# Populate matrix with count of number of genes with padj <=0.05 and log2fc >2 for increase or <-2 for decrease
for (i in temp.DGE.comparisons) {
  temp <- as.data.frame(results(dds, name=i)) %>% tibble::rownames_to_column("gene_id") # best df <- matrix
  temp.increase <-  temp[temp$padj <=0.05 & temp$log2FoldChange > 2 & !is.na(temp$padj), ]
  temp.decrease <-  temp[temp$padj <=0.05 & temp$log2FoldChange < -2 & !is.na(temp$padj), ]
  temp.DGE.count[i, "increase"] <- nrow(temp.increase)
  temp.DGE.count[i, "decrease"] <- nrow(temp.decrease)
}
rm(i)
print(temp.DGE.count)

# Export Table after removing Intercept
temp.DGE.count <- temp.DGE.count[!grepl('Intercept', rownames(temp.DGE.count)),, drop = FALSE ] 
write.csv(temp.DGE.count, file = paste0('PCA/DEG_count_', timepoint, '.csv'))

# Filter matrix for timepoint / -Intercept
# alter names with gsub
# filters on matching timepoint comparison string, ie treatment_RadioTherapy_{14d_vs_vehicle_14d}	
includes=paste0(timepoint,'_vs_',control_name,'_',timepoint)
temp.DGE.count2 <- temp.DGE.count[grepl(includes, rownames(temp.DGE.count)),, drop = FALSE ] # MUST ADD ,,DROP=FALSE WHEN SINGLE RESULT OTHERWISE DROPS ROWNAMES
print(temp.DGE.count2)
rm(includes)
rownames(temp.DGE.count2) <- gsub(pattern = 'treatment_', replacement = '', x = rownames(temp.DGE.count2))

# DF, gather, plot 
temp.DGE.count2 <- as.data.frame(temp.DGE.count2) %>% tibble::rownames_to_column("Comparison") # best df <- matrix
temp.DGE.count3 <- temp.DGE.count2 %>% gather('increase', 'decrease', key = "change", value = "DEG_count")

ggplot(temp.DGE.count3, aes(fill=change, y=DEG_count, x=Comparison)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title=paste0('Differentially Expressed Genes', ' ', timepoint),
       subtitle = 'cut off: adj-pvalue 0.05, Log2FC +/- 2',
       x=element_blank()) +
  scale_fill_manual(values = c("increase" = "#DC000099", "decrease" = "#3C548899")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300), breaks = seq(0, 300, by = 100)) +
  theme(plot.title = element_text(color = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
        plot.subtitle = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color = "black", size = 10, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        rect = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5, linetype = "solid", fill = NA)
        ) 
       
pdf(paste0('PCA/DEG_count_', timepoint, '.pdf'), width=5, height=7)
print(last_plot())
dev.off()

#cleanup
rm(list = ls(pattern = 'temp'))
rm(dds, timepoint)

}
