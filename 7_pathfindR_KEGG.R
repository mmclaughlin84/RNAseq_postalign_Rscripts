######################################################################
# __________ pathfindR  Pathway Analysis__________ #
######################################################################

# https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html
# https://cran.r-project.org/web/packages/pathfindR/vignettes/visualization_vignette.html
# https://cran.r-project.org/web/packages/pathfindR/vignettes/non_hs_analysis.html

rm(list = ls()) #clear environment
setwd("/Users/mclaughinm/Desktop/RNAseq_EPMOC2/analysis") # set working directory as *ANALYSIS* folder


dev.off()

getwd() # get working directory
list.files()
dir.create('pathfindR')
setwd("/Users/mclaughinm/Desktop/RNAseq_EPMOC2/analysis") # set working directory as *ANALYSIS* folder

library(pathfindR)
library(tidyverse)
library(org.Hs.eg.db)
library(KEGGREST)
library(KEGGgraph)
library(biomaRt)


######################################################################
# 1) CREATE FILE LIST OF DEGs
######################################################################

# List timepoints
# Create blank matrix with number of comparisons as nrows needed
timepoints=c('D3', 'D10')
pathfindr_file_list <- data.frame(matrix(NA, nrow = 3, ncol = 0))

for (timepoint in timepoints) {
pathfindr_file_path=paste0('DESeq_export/', timepoint) # Indicate folder by timepoint
temp_file <- data.frame(file = c(list.files(path=paste0(pathfindr_file_path, '/')))) # list files in folder
temp_file <- temp_file[grepl(paste0(timepoint,'_vs_'), temp_file$file), , drop = FALSE] # restrict to timepoint
pathfindr_file_list$file <- temp_file$file # move column to created dataframe now that the length is correct
colnames(pathfindr_file_list)[colnames(pathfindr_file_list) == 'file'] <- paste0(timepoint) # rename by timepoint
print(pathfindr_file_list)
rm(temp_file, timepoint, pathfindr_file_path)
}


######################################################################
# 2) pathfindR ANALYSIS
######################################################################

for (timepoint in timepoints) {
  
  # Extract list of files by timepoint to pass to nested second loop
  print(file_list <- pathfindr_file_list[[timepoint]])
  
  for (file in file_list) {
    
    # IMPORT DEG.csv
    # Filter on padj and log fold change (remove NA needed) and extract only 
    # 3 Column in put required. MGI, log2FC, adjpvalue
    pathfindr_imported_DEGs <- read.csv(paste0('DESeq_export/', timepoint, '/', file)) %>% dplyr::rename(Gene.symbol = 1) # this will not work if read_csv instead of read.csv
    pathfindr_imported_DEGs <- pathfindr_imported_DEGs[((pathfindr_imported_DEGs$log2FoldChange >= +2 | pathfindr_imported_DEGs$log2FoldChange <= -2) &pathfindr_imported_DEGs$padj <= 0.05 & pathfindr_imported_DEGs$padj > 0 & !is.na(pathfindr_imported_DEGs$padj)), c('Gene.symbol', 'log2FoldChange', 'padj')]
    colnames(pathfindr_imported_DEGs) <- c('Gene.symbol', 'logFC', 'adj.P.Val')
    head(pathfindr_imported_DEGs[order(pathfindr_imported_DEGs$adj.P.Val),],10)
    
    ###### RUNNING PATHFINDR
    # Protein-Protein Interaction Networks (PIN)
    # Default = Biogrid, but KEGG, GeneMania and IntAct also possible
    # 1) PIN (default - biogrid), 2) Algorithm (default - greedy), 3) Geneset (default - KEGG)
    ?run_pathfindR
    pathfindr_output <- run_pathfindR(input = pathfindr_imported_DEGs,
                                      gene_sets = "mmu_KEGG", # Ms gene set = mmu_KEGG
                                      pin_name_path = "mmu_STRING", # Ms PIN = mmu_STRING
                                      search_method = "GR", # GR=greedy, SA=annealling, GA=genetic [3 Algorithms Greedy, Annealing, Genetic]
                                      convert2alias = FALSE, # only TRUE for Hm
                                      p_val_threshold = 0.05,
                                      min_gset_size = 10, # 10 = default
                                      max_gset_size = 300, # 300 = default
                                      adj_method = "bonferroni", # default = bonferroni not fdr
                                      enrichment_threshold = 0.05, # default = 0.05
                                      output_dir = paste0('pathfindR/', file)
    )
    assign(paste0('pathfindR_result_', file), pathfindr_output)
    rm(pathfindr_imported_DEGs, pathfindr_output, file)
    
  }
  rm(timepoint, file_list)
}

#save.image(file = "pathfindR/environment.RData")
#load('pathfindR/environment.Rdata')


######################################################################
# 3) EXPORTING INDIVIDUAL PathfindR PLOTS
######################################################################

# List all outputs from pathfindR in steps above
(pathfindR_output_list <- ls(pattern = 'pathfindR_result_'))

for (output in pathfindR_output_list) {

# manual
output = 'pathfindR_result_treatment_RT_ATRi_D3_vs_Vehicle_D3.csv'
manual_title='ATRi/RT Day 3'
##### For Test Purposes
##### (output <- pathfindR_output_list[4])

# Abbreviate output name to comparison name for title of plot and export file name
output_abbreviated <- gsub(pattern = 'pathfindR_result_treatment_', replacement = '', output)
output_abbreviated <- gsub(pattern = '.csv', replacement = '', output_abbreviated)

# Enrichment chart is a ggplot object
plot <- enrichment_chart(result_df = get(output), top_terms = 10) +
  labs(title = manual_title,
       x = "Fold Enrichment", #expression to insert mathematical symbols // Use ~ in place of space, will auto convert greek
       y = NULL) + # x = NULL or y = NULL removes axis label
  #scale_x_continuous(expand = c(0,0), limits = c(0, 10), breaks = seq(0, 10, by = 2), minor_breaks = NULL) +
  theme(aspect.ratio=4/1.5) +
  theme(plot.title = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 1, face = "bold", family = "sans"),
        axis.text.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.y = element_text(color = "black", size = 7, angle = 90, vjust = 0.5, hjust = 0.5, face = "plain", family = "sans")
        ) +
  #theme(legend.key.size = unit(.4, "cm")) +
  guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3, title.theme = element_text(size = 7), label.theme = element_text(size = 7))) +
  guides(size = guide_legend(keyheight = 0.8, title = "#genes", legend.key.size = 10, title.theme = element_text(size = 7), label.theme = element_text(size = 7))) +
  scale_size(range = c(1, 4),)  # adjusts scale of bubble

print(plot)

pdf(paste0('pathfindr/', 'fig_', output_abbreviated, '.pdf'), width=3.5, height=2.1)
print(plot)
dev.off()

}


######################################################################
# 4) CLUSTERED PathfindR PLOTS
######################################################################

### Skipping Clustered pathfindR plots as not usful at this time
#RA_clustered <- cluster_enriched_terms(pathfindr.output, plot_dend = FALSE, plot_clusters_graph = FALSE)
#enrichment_chart(RA_clustered, plot_by_cluster = TRUE)
#RA_selected <- subset(RA_clustered, Cluster %in% 1:8)
#enrichment_chart(RA_selected, plot_by_cluster = TRUE)

######################################################################
# 5) GROUPED PathfindR PLOTS
######################################################################

### Grouped PathfindR plots also not useful
# List all outputs from pathfindR in steps above
#(pathfindR_output_list <- ls(pattern = 'pathfindR_result_'))

# Combining Results
#combined_df <- combine_pathfindR_results(result_A = pathfindR_result_treatment_ATRi_D3_vs_Vehicle_D3.csv, 
#                                         result_B = pathfindR_result_treatment_RT_D3_vs_Vehicle_D3.csv,
#                                         plot_common = FALSE)
#combined_results_graph(combined_df,
#                       selected_terms = "common",
#                       use_description = FALSE,
#                       layout = "stress",
#                       node_size = "num_genes")

# You may run `combined_results_graph()` to create visualizations of combined term-gene graphs of selected terms

######################################################################
# 6a) FIGURE ASSEMBLY DAY 3
######################################################################

ls(pattern = 'pathfindR_result_treatment')
rm(list = ls(pattern = 'D10_vs_|temp|mMCP|go_result|heatmap|topgo',), 'h', 'plot')
ls(pattern = 'pathfindR_result_treatment')


# Column names
# Copy only needed columns to dataframe to be tidied
# Add blank count column
ls(pattern = 'pathfindR_result')
colnames(pathfindR_result_treatment_ATRi_D3_vs_Vehicle_D3.csv)
tidy_ATRi_D3 <- pathfindR_result_treatment_ATRi_D3_vs_Vehicle_D3.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_RT_D3 <- pathfindR_result_treatment_RT_D3_vs_Vehicle_D3.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_ATRiRT_D3 <- pathfindR_result_treatment_RT_ATRi_D3_vs_Vehicle_D3.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_ATRi_D3$count <- NA
tidy_RT_D3$count <- NA
tidy_ATRiRT_D3$count <- NA

# tidy_ATRi_D3
for (rownumber in (1:nrow(tidy_ATRi_D3))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRi_D3[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRi_D3[rownumber,'Down_regulated']), ","))))
  tidy_ATRi_D3[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}
# tidy_RT_D3
for (rownumber in (1:nrow(tidy_RT_D3))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_RT_D3[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_RT_D3[rownumber,'Down_regulated']), ","))))
  tidy_RT_D3[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}
# tidy_ATRiRT_D3
for (rownumber in (1:nrow(tidy_ATRiRT_D3))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRiRT_D3[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRiRT_D3[rownumber,'Down_regulated']), ","))))
  tidy_ATRiRT_D3[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}

# Select only pathways with at least 4 genes present
tidy_ATRi_D3 <- tidy_ATRi_D3[tidy_ATRi_D3$count >= 4,]
tidy_RT_D3 <- tidy_RT_D3[tidy_RT_D3$count >= 4,]
tidy_ATRiRT_D3 <- tidy_ATRiRT_D3[tidy_ATRiRT_D3$count >= 4,]

# Export For Supplementary
dir.create('pathfindR/sup_table/')
write.csv(tidy_ATRi_D3, file = 'pathfindR/sup_table/ATRi_D3.csv')
write.csv(tidy_RT_D3, file = 'pathfindR/sup_table/RT_D3.csv')
write.csv(tidy_ATRiRT_D3, file = 'pathfindR/sup_table/ATRiRT_D3.csv')


# Select only fold enrichment and pvalue 
colnames(tidy_ATRi_D3)
tidy_ATRi_D3_merge <- tidy_ATRi_D3[,c('Term_Description','Fold_Enrichment','lowest_p')]
tidy_RT_D3_merge <- tidy_RT_D3[,c('Term_Description','Fold_Enrichment','lowest_p')]
tidy_ATRiRT_D3_merge <- tidy_ATRiRT_D3[,c('Term_Description','Fold_Enrichment','lowest_p')]

# Add treatment column
tidy_ATRi_D3_merge$treatment <- 'ATRi'
tidy_RT_D3_merge$treatment <- 'RT'
tidy_ATRiRT_D3_merge$treatment <- 'ATRi/RT'

# Merge all as tall data with treatment column
tidy_all_merged_d3 <- merge(tidy_ATRi_D3_merge, tidy_RT_D3_merge, all = TRUE) %>% merge(., tidy_ATRiRT_D3_merge, all = TRUE)

### *** CORRECTED TOP 25 TERMS BY P VALUE *** ###
tidy_all_merged_selected <- tidy_all_merged_d3[order(tidy_all_merged_d3$lowest_p),] # sort by pvalue
tidy_all_merged_selected$Term_Description <- factor(tidy_all_merged_selected$Term_Description, levels = unique(tidy_all_merged_selected$Term_Description)) # terms to factors
tidy_all_merged_selected$treatment <- factor(tidy_all_merged_selected$treatment, levels = c('ATRi', 'RT','ATRi/RT')) # treatment to factors
tidy_all_merged_selected <- tidy_all_merged_selected[tidy_all_merged_selected$Term_Description %in% head(levels(tidy_all_merged_selected$Term_Description), 15),] # top 18 terms to match old graph
str(tidy_all_merged_selected)

# 'Bottom' 25 ***OLD - DELETE***
#tidy_all_merged_selected_d3 <- head(tidy_all_merged_d3[order(tidy_all_merged_d3$lowest_p),],25)
# Convert 1) Terms and 2) Treatments to factors
#tidy_all_merged_selected_d3$Term_Description <- factor(tidy_all_merged_selected_d3$Term_Description, levels = unique(tidy_all_merged_selected_d3$Term_Description))
#tidy_all_merged_selected_d3$treatment <- factor(tidy_all_merged_selected_d3$treatment, levels = c('ATRi', 'RT','ATRi/RT'))
#str(tidy_all_merged_selected_d3$treatment)


# ggplot
colnames(tidy_all_merged_selected)
plot_merged <- ggplot(tidy_all_merged_selected, aes(x=log10(lowest_p), y=reorder(Term_Description, desc(Term_Description)), size=Fold_Enrichment, color=treatment)) +
  geom_point(alpha=1) +
  scale_color_manual(values = c('#b855d8','#08bf01','#5657f9')) +
  scale_x_reverse() +
  scale_size_binned(range = c(0, 5)) + # come back
  facet_wrap(vars(treatment), nrow = 1, ncol = 3, ) +
  labs(title = 'Pathway Analysis (pathfindR) Day 3',
       x = "Log10 P-value", #expression to insert mathematical symbols // Use ~ in place of space, will auto convert greek
       y = NULL) + # x = NULL or y = NULL removes axis label
  theme(plot.title = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = -4, face = "bold", family = "sans"),
        axis.text.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.text.y = element_text(color = "black", size = 7, angle = 0, hjust = 1, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.y = element_text(color = "black", size = 7, angle = 90, vjust = 0.5, hjust = 0.5, face = "plain", family = "sans"),
        strip.text = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = -1, face = "plain", family = "sans"),
        strip.background = element_rect(fill='white', colour=NULL,size=1),
        panel.spacing = unit(1, "mm"),
        legend.spacing = unit(2, "mm"),
        legend.margin = margin(c(0,3,0,-6)),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill='white'),
        panel.border = element_rect(fill = NA, colour = 'black', size = 1),
        panel.grid.major = element_line(colour = '#e9eef0', size = 0.25)
        ) +
  guides(colour = guide_legend(order = 1, keyheight = 0.8, title = "Treatment", legend.key.size = 8, title.theme = element_text(size = 7), label.theme = element_text(size = 7))) +
  guides(size = guide_legend(order = 2, keyheight = 0.8, title = "Fold Enrichment", legend.key.size = 8, title.theme = element_text(size = 7), label.theme = element_text(size = 7)))

print(plot_merged)


pdf(paste0('pathfindr/', 'updated_merged_plot_', 'Day3', '.pdf'), width=5, height=2.5)
print(plot_merged)
dev.off()



######################################################################
# 6b) FIGURE ASSEMBLY DAY 10
######################################################################

ls(pattern = 'pathfindR_result_treatment')
rm(list = ls(pattern = 'D3_vs_|temp|mMCP|go_result|heatmap|topgo',), 'h', 'plot')
ls(pattern = 'pathfindR_result_treatment')


# Column names
# Copy only needed columns to dataframe to be tidied
# Add blank count column
colnames(pathfindR_result_treatment_ATRi_D10_vs_Vehicle_D10.csv)
tidy_ATRi_D10 <- pathfindR_result_treatment_ATRi_D10_vs_Vehicle_D10.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_RT_D10 <- pathfindR_result_treatment_RT_D10_vs_Vehicle_D10.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_ATRiRT_D10 <- pathfindR_result_treatment_RT_ATRi_D10_vs_Vehicle_D10.csv[,c('Term_Description','Fold_Enrichment','lowest_p','Up_regulated','Down_regulated')]
tidy_ATRi_D10$count <- NA
tidy_RT_D10$count <- NA
tidy_ATRiRT_D10$count <- NA

# tidy_ATRi_D10
for (rownumber in (1:nrow(tidy_ATRi_D10))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRi_D10[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRi_D10[rownumber,'Down_regulated']), ","))))
  tidy_ATRi_D10[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}
# tidy_RT_D10
for (rownumber in (1:nrow(tidy_RT_D10))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_RT_D10[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_RT_D10[rownumber,'Down_regulated']), ","))))
  tidy_RT_D10[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}
# tidy_ATRiRT_D10
for (rownumber in (1:nrow(tidy_ATRiRT_D10))) {
  count1 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRiRT_D10[rownumber,'Up_regulated']), ","))))
  count2 <- data.frame(table(unlist(strsplit(as.character(tidy_ATRiRT_D10[rownumber,'Down_regulated']), ","))))
  tidy_ATRiRT_D10[rownumber,'count'] <- (sum(count1$Freq) + sum(count2$Freq))
}

# Select only pathways with at least 4 genes present
tidy_ATRi_D10 <- tidy_ATRi_D10[tidy_ATRi_D10$count >= 4,]
tidy_RT_D10 <- tidy_RT_D10[tidy_RT_D10$count >= 4,]
tidy_ATRiRT_D10 <- tidy_ATRiRT_D10[tidy_ATRiRT_D10$count >= 4,]

# Export For Supplementary
dir.create('pathfindR/sup_table/')
write.csv(tidy_ATRi_D10, file = 'pathfindR/sup_table/ATRi_D10.csv')
write.csv(tidy_RT_D10, file = 'pathfindR/sup_table/RT_D10.csv')
write.csv(tidy_ATRiRT_D10, file = 'pathfindR/sup_table/ATRiRT_D10.csv')

# Select only fold enrichment and pvalue 
colnames(tidy_ATRi_D10)
tidy_ATRi_D10_merge <- tidy_ATRi_D10[,c('Term_Description','Fold_Enrichment','lowest_p')]
tidy_RT_D10_merge <- tidy_RT_D10[,c('Term_Description','Fold_Enrichment','lowest_p')]
tidy_ATRiRT_D10_merge <- tidy_ATRiRT_D10[,c('Term_Description','Fold_Enrichment','lowest_p')]

# Add treatment column
tidy_ATRi_D10_merge$treatment <- 'ATRi'
tidy_RT_D10_merge$treatment <- 'RT'
tidy_ATRiRT_D10_merge$treatment <- 'ATRi/RT'

# Merge all as tall data with treatment column
tidy_all_merged <- merge(tidy_ATRi_D10_merge, tidy_RT_D10_merge, all = TRUE) %>% merge(., tidy_ATRiRT_D10_merge, all = TRUE)

### *** CORRECTED TOP 25 TERMS BY P VALUE *** ###
tidy_all_merged_selected <- tidy_all_merged[order(tidy_all_merged$lowest_p),] # sort by pvalue
tidy_all_merged_selected$Term_Description <- factor(tidy_all_merged_selected$Term_Description, levels = unique(tidy_all_merged_selected$Term_Description)) # terms to factors
tidy_all_merged_selected$treatment <- factor(tidy_all_merged_selected$treatment, levels = c('ATRi', 'RT','ATRi/RT')) # treatment to factors
tidy_all_merged_selected <- tidy_all_merged_selected[tidy_all_merged_selected$Term_Description %in% head(levels(tidy_all_merged_selected$Term_Description), 15),] # top 18 terms to match old graph
str(tidy_all_merged_selected)

# 'Bottom' 25 ***OLD - DELETE LATER***
#tidy_all_merged_selected <- head(tidy_all_merged[order(tidy_all_merged$lowest_p),],25)
# Convert 1) Terms and 2) Treatments to factors
#tidy_all_merged_selected$Term_Description <- factor(tidy_all_merged_selected$Term_Description, levels = unique(tidy_all_merged_selected$Term_Description))
#tidy_all_merged_selected$treatment <- factor(tidy_all_merged_selected$treatment, levels = c('ATRi', 'RT','ATRi/RT'))
#str(tidy_all_merged_selected$treatment)


# ggplot
colnames(tidy_all_merged_selected)
plot_merged <- ggplot(tidy_all_merged_selected, aes(x=log10(lowest_p), y=reorder(Term_Description, desc(Term_Description)), size=Fold_Enrichment, color=treatment)) +
  geom_point(alpha=1) +
  scale_color_manual(values = c('#b855d8','#08bf01','#5657f9')) +
  scale_x_reverse() +
  scale_size_binned(range = c(0, 3), n.breaks = 4) + # come back
  facet_wrap(vars(treatment), nrow = 1, ncol = 3, ) +
  labs(title = 'Pathway Analysis (pathfindR) Day 10',
       x = "Log10 P-value", #expression to insert mathematical symbols // Use ~ in place of space, will auto convert greek
       y = NULL) + # x = NULL or y = NULL removes axis label
  theme(plot.title = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = -4, face = "bold", family = "sans"),
        axis.text.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.x = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.text.y = element_text(color = "black", size = 7, angle = 0, hjust = 1, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.y = element_text(color = "black", size = 7, angle = 90, vjust = 0.5, hjust = 0.5, face = "plain", family = "sans"),
        strip.text = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = -1, face = "plain", family = "sans"),
        strip.background = element_rect(fill='white', colour=NULL,size=1),
        panel.spacing = unit(1, "mm"),
        legend.spacing = unit(2, "mm"),
        legend.margin = margin(c(0,3,0,-6)),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill='white'),
        panel.border = element_rect(fill = NA, colour = 'black', size = 1),
        panel.grid.major = element_line(colour = '#e9eef0', size = 0.25)
  ) +
  guides(colour = guide_legend(order = 1, keyheight = 0.8, title = "Treatment", legend.key.size = 8, title.theme = element_text(size = 7), label.theme = element_text(size = 7))) +
  guides(size = guide_legend(order = 2, keyheight = 0.8, title = "Fold Enrichment", legend.key.size = 8, title.theme = element_text(size = 7), label.theme = element_text(size = 7)))

print(plot_merged)


pdf(paste0('pathfindr/', 'updated_merged_plot_', 'day10', '.pdf'), width=5.8, height=2.5)
print(plot_merged)
dev.off()





