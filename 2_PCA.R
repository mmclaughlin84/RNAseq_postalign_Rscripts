#rm(list = ls())

################################################ 
#_____Basic QC: PCA/Heatmap/#DEGs_____ #
################################################

### Run to end of DEseq analysis and switch to this script

# 1) PCA (used normalised genes - timepoint agnostic)
# 2) Heatmap All Genes (used normalised genes - timepoint agnostic)
# 3) DEG Count (!!! requires timepoint and control_name details from DESeq !!!)


############################
# 1) PCA
############################

### !!!NOTE!!! 3 lines will need edited for treatment/timepoint variables if different from
# $timepoint
# $treatment_wo_time

# This has been switched to import normalised counts from csv file to avoid that dds files do not have a consistent naming structure
pcaCounts <- data.frame(read.csv('DESeq_export/normalised_counts.csv', row.names = 1)) %>% tibble::rownames_to_column("gene_id")

# PCA: Filter Normalised Counts from DESeq dds dataset (rows with zero variance cause error)
# based on https://www.youtube.com/watch?v=0Jp4gsfOLMs&t=364s
rownames(pcaCounts) <- pcaCounts$gene_id # add rownames from gene id column
pcaCounts <- pcaCounts[ ,-1] # delete gene_id column
pcaCounts <- as.matrix(pcaCounts) # make matrix
keep <- rowSums(pcaCounts) >= 340 # results in a TRUE/FALSE output # set at 10 reads per sample. Has big impact!
pcaCounts <- pcaCounts[keep,] # only keeps rows = TRUE
rm(keep)

# PCA: prcomp expects samples as rows and genes as columns
# PCA: plot pc1 and pc2
pca <- prcomp(t(pcaCounts), scale=TRUE)
summary(pca)
plot(pca$x[,1], pca$x[,2]) # basic plot
pca.plot.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1],Y=pca$x[,2])

# add timepoint and treatment from colData
pca.plot.data$timepoint <- colData$timepoint
pca.plot.data$treatment_wo_time <- colData$treatment_wo_time
pca.plot.data

# PCA: scree plot, how much variation is explained by each principle component
pca.var <- pca$sdev^2 # sq of std dev of each pc = variance
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # calc as %
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")


# PCA Plot by ggplot
ggplot(data=pca.plot.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(color = treatment_wo_time, shape = timepoint)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")

dir.create("PCA/", showWarnings = FALSE)
ggsave('PCA/PCA_plot.pdf', plot = last_plot(), width = 15, height = 10, units = 'cm', dpi = 300)


## PCA: get the top 10 genes that contribute most to pc1
pca.loading_scores <- pca$rotation[,1]
pca.gene_scores <- abs(pca.loading_scores) ## get the magnitudes
pca.gene_score_ranked <- sort(pca.gene_scores, decreasing=TRUE)
(pca.top_10_genes <- names(pca.gene_score_ranked[1:10])) ## show the names of the top 10 genes
pca$rotation[pca.top_10_genes,1] ## show the scores (and +/- sign)

# cleanup
normalisedCounts <- pcaCounts # 1) matrix 2) filtered - for use in heatmap below
rm(list = ls(pattern = 'pca'))


############################
# 2) Heatmap All Genes
############################

##### 6c) QC: Heatmap (All Normalised Counts) #####

# Reuse pcaCounts matrix from PCA prcomp step
normalisedCounts.heatmap <- t(normalisedCounts)
normalisedCounts.heatmap <- scale(normalisedCounts.heatmap)

heatmap_plot <- Heatmap(t(normalisedCounts.heatmap),
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
pdf('PCA/heatmap_all.pdf', width=5, height=5)
draw(heatmap_plot)
dev.off()

#cleanup
rm(list = c('heatmap_plot', 'normalisedCounts', 'normalisedCounts.heatmap'))


############################
# 3) DEGs Count
############################

# run timepoint and control_name lines in DESeq_timepoints file (lines 24-25ish)

for (timepoint in timepoints) {

# Time timepoint specific dds output and convert it to dds naming to ease loop script
dds_by_timepoint=paste0('dds_', timepoint)
dds='dds'
assign(dds, get(dds_by_timepoint))
rm(dds_by_timepoint)

  
# Get comparison list from resultsNames(dds)
# Create matrix 
(temp.DGE.comparisons=resultsNames(dds)) # lookup results for next line
temp.DGE.count <- matrix(nrow = length(temp.DGE.comparisons), ncol = 2)
rownames(temp.DGE.count) <- temp.DGE.comparisons
colnames(temp.DGE.count) <- c('increase', 'decrease')

# Populate Matrix
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1400), breaks = seq(0, 2000, by = 200)) +
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
