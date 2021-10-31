
### RUN DESEQ2 up to export of files with adjpvalues ###
### Script setup to import these files ###


#_________________GSEA___________________#
### After DEseq (ensure generic input) ###

library(fgsea)
library(reshape2)  # needed?
library(circlize) # needed?

# File List
# This is spectacularly hamfisted
GSEA_file_path='DEseq_export_14d' # Populate folder and it should do the rest
GSEA_file_list <- data.frame(file=c(list.files(path=paste0(GSEA_file_path, '/')))) 
GSEA_file_list <- GSEA_file_list[grepl('_vs_', GSEA_file_list$file), ] 
GSEA_file_list <- data.frame(GSEA_file_list)
names(GSEA_file_list)[1] <- 'file'
GSEA_file_list

# Manual setup of for loop 
i='treatment_RadioTherapy_BI880_21d_vs_vehicle_21d.csv'

for (i in GSEA_file_list$file) {

# [GSEA] 1) Import DEG.csv
# gene_id must be HNGC(Hm) or MGI(Mm) and all caps
gsea.deg.csv <- read_csv(paste0(GSEA_file_path, '/', i)) %>% dplyr::rename(gene_id = 1)
gsea.deg.csv <- gsea.deg.csv[(gsea.deg.csv$log2FoldChange > 2 | gsea.deg.csv$log2FoldChange < -2) & gsea.deg.csv$padj <= 0.05 & !is.na(gsea.deg.csv$padj), ] 
gsea.deg.csv$gene_id <- toupper(gsea.deg.csv$gene_id)

# remove the NAs, averaging statistics for a multi-hit symbol
# rank genes by DEseq stat variable
gsea.deg.csv2 <- gsea.deg.csv %>% dplyr::select(gene_id, stat) %>% na.omit() %>% distinct() %>% group_by(gene_id) %>% summarize(stat=mean(stat))
ranks <- gsea.deg.csv2$stat
names(ranks) <- gsea.deg.csv2$gene_id

# Load gene set (all hallmarks) https://www.gsea-msigdb.org/gsea/index.jsp
pathways.hallmark <- gmtPathways("GSEA/h.all.v7.2.symbols.gmt")

#Running fgsea algorithm: # Tidy  results:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) # order by normalized enrichment score (NES)

# Maps hallmark pathways to gene_ids
gene.in.pathway <- pathways.hallmark %>% enframe("pathway", "gene_id") %>% unnest(cols = c(gene_id)) %>% inner_join(gsea.deg.csv, by="gene_id")


# Plot  normalized enrichment scores
# ifelse indicating significant as red
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
gsea.bar.plot <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = cols) +
  theme_bw() +
 theme(plot.title = element_text(color = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.text.x = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.x = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain", family = "sans"),
        axis.text.y = element_text(color = "black", size = 7, angle = 0, hjust = 1, vjust = 0.5, face = "plain", family = "sans"),
        axis.title.y = element_text(color = "black", size = 7, angle = 0, vjust = 0.5, hjust = 0.5, face = "plain", family = "sans"),
        axis.line = element_line(colour = "black", size = .25, linetype = "solid")) +
  labs(x=element_blank(),
       y="Normalized Enrichment Score",
       title=paste0(i,' GSEA'))

pdf(paste0('GSEA/', i, '.pdf'), width=9, height=7)
print(gsea.bar.plot)
dev.off()

}





#__________ Specific Enrichment  Plot_______#
# Enrichment plot for E2F target gene set
plotEnrichment(pathway = pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]], ranks, ticksSize = 0.2) +
  labs(title='HALLMARK_INFLAMMATORY_RESPONSE')




# Removal of Excessive Environment Objects
rm(ranks)
rm(pathways.hallmark)
rm(fgseaRes)
rm(fgseaResTidy)
rm(gene.in.pathway)
rm(cols)













##### DID NOT GO BEYOND THIS POINT #####



plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
              gseaParam=0.5)

#________ Heatmap Plot_____________#
# pathways with significant enrichment score
sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
sig.gen <- unique(na.omit(gene.in.pathway$gene_id[gene.in.pathway$pathway %in% sig.path]))

### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
h.dat <- dcast(gene.in.pathway[, c(1,2)], gene_id~pathway)
rownames(h.dat) <- h.dat$gene_id
h.dat <- h.dat[, -1]

h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]

# keep those genes with 3  or more occurnes
table(data.frame(rowSums(h.dat)))

# 1       2    3    4    5    6 
# 1604  282   65   11    1    1 
h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]

#
topTable <- res[res$gene_id %in% rownames(h.dat), ]
rownames(topTable) <- topTable$gene_id
# match the order of rownames in toptable with that of h.dat
topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(h.dat))

# colour bar for -log10(adjusted p-value) for sig.genes
dfMinusLog10FDRGenes <- data.frame(-log10(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'padj']))
dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0

# colour bar for fold changes for sigGenes
dfFoldChangeGenes <- data.frame(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'log2FoldChange'])

# merge both
dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                         ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
colours <- list(
  'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
haGenes <- rowAnnotation(
  df = dfGeneAnno,
  col = colours,
  width = unit(1,'cm'),
  annotation_name_side = 'top')

# Now a separate color bar for the GSEA enrichment padj. This will 
# also contain the enriched term names via annot_text()

# colour bar for enrichment score from fgsea results
dfEnrichment <- fgseaRes[, c("pathway", "NES")]
dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
dd <- dfEnrichment$pathway
dfEnrichment <- dfEnrichment[, -1]
rownames(dfEnrichment) <- dd
colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
haTerms <- HeatmapAnnotation(
  df = dfEnrichment,
  Term = anno_text(
    colnames(h.dat),
    rot = 45,
    just = 'right',
    gp = gpar(fontsize = 12)),
  annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
  annotation_name_side = 'left')
# now generate the heatmap
hmapGSEA <- Heatmap(h.dat,
                    name = 'GSEA hallmark pathways enrichment',
                    split = dfGeneAnno[,2],
                    col = c('0' = 'white', '1' = 'forestgreen'),
                    rect_gp = gpar(col = 'grey85'),
                    cluster_rows = TRUE,
                    show_row_dend = TRUE,
                    row_title = 'Top Genes',
                    row_title_side = 'left',
                    row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_title_rot = 90,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_names_side = 'left',
                    row_dend_width = unit(35, 'mm'),
                    cluster_columns = TRUE,
                    show_column_dend = TRUE,
                    column_title = 'Enriched terms',
                    column_title_side = 'top',
                    column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                    column_title_rot = 0,
                    show_column_names = FALSE,
                    show_heatmap_legend = FALSE,
                    clustering_distance_columns = 'euclidean',
                    clustering_method_columns = 'ward.D2',
                    clustering_distance_rows = 'euclidean',
                    clustering_method_rows = 'ward.D2',
                    bottom_annotation = haTerms)

tiff("GSEA_enrichment_2.tiff", units="in", width=13, height=22, res=400)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

