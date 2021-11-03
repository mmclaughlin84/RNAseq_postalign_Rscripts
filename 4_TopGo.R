###########################
### GO Enrichment w/ TopGo
###########################

# Set Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # sets working directory based on script location - should be in analysis folder
getwd() # check

# Load libraries / create export directory
library(GO.db)
library(topGO)
library(biomaRt)
library(Rgraphviz) # needed for (not very useful) GO plots
library(ComplexHeatmap)
library(DESeq2)
suppressWarnings(dir.create('TopGO'))

# TopGo Overviews
# https://bioconductor.org/packages/release/bioc/html/topGO.html
# https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/

# Load saved DESeq2 script environment
load(file = "DESeq_export/DESeq_output.RData") # Avoids rerunning DESeq2 - is >50Mb so if offsite on VPN can take a moment to load

# Note: There is a line at the start of the second for-loop defining a 'Gene of interest' as a DEG <2 or >2 log2fc plus <= 0.05 padj

#########################################################
# STEP 1: Import 'Gene Universe' and Compile GO terms
#########################################################

# Gene Universe = all genes listed in the normalised_counts.csv file
geneUniverse <- data.frame(read.csv('DESeq_export/normalised_counts.csv', row.names = 1)) %>% tibble::rownames_to_column("gene_id")
geneUniverse <- as.character(geneUniverse[,'gene_id']) # wrap in as.character() as a precaution
print(head(geneUniverse))

# Using BioMart look up GO terms for all genes in geneUniverse - lookup takes a while
# go_id = GO term accession
# namespace_1003 = GO domain (biological_process/molecular_function/cellular_component)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #  takes a moment - often says 'Ensembl site unresponsive, trying ...' which is normal
head(listAttributes(ensembl), 5) # should you want to check the filter/attributes categories in the database
go_ids = getBM(attributes = c("go_id", "mgi_symbol", "namespace_1003"), # namespace_1003 (GO domain) in place of name_1006 (GO term name)
               filters = 'mgi_symbol', 
               values = geneUniverse, 
               mart = ensembl)
rm(ensembl)

# converts go_id list to wide unstacked format with go_terms in single row with gene_id. Gene ids with no GO terms, they appear as ""
# Later loaded as 'gene_2_GO' in TopGo object
go_ids_wide <- unstack(go_ids[,c(1,2)])


#########################################################
# Step 2: Import Genes of Interest and Make TopGo Object
#########################################################
# FOR-LOOP OF DESeq2 "TIMEPOINTS": START
for (go_timepoint in timepoints) {
  
  # [2a Create Timepoint linked dds]: Assign timepoint specific dds to 'dds'
  cat(magenta(paste0('Loading Timepoint Specific DESeq2 DDS: ', paste0('dds_', go_timepoint,'\n'))))
  dds <- get(paste0('dds_', go_timepoint))
  
  # [2b Results comparisons lookup]: Need to pull from timepoint specific dds object filter with grepl so only within timepoint comparisons
  cat(magenta(paste0('Extract DESeq2 DDS Comparisons List:','\n'))) # Just printing to console
  go_timepoint_comparisons <- resultsNames(dds) # unfiltered DESeq2 comparisons list
  print(data.frame(unfiltered_comparisons = go_timepoint_comparisons)) # Just printing to console
  keep <- grepl(paste0(go_timepoint, '_vs_' ), go_timepoint_comparisons) # TRUE/FALSE Filter with 'timepoint'_vs_ as pattern
  go_timepoint_comparisons <- go_timepoint_comparisons[keep]
  cat(magenta(paste0('Filtered Comparisons List:','\n'))) # Just printing to console
  print(data.frame(FILTERED_comparisons = go_timepoint_comparisons)) # Just printing to console
  # NOTE: Groups without any deferentially expressed genes will result in an error - need to be removed with the below silenced code
  #keep <- !grepl(paste0('treatment_GROUPTOREMOVE_', go_timepoint), go_timepoint_comparisons)
  #(go_timepoint_comparisons <- go_timepoint_comparisons[keep])
  rm(keep)
  
  #!!! FOR-LOOP OF DESeq2 "COMPARISONS" WITHIN EACH TIMEPOINT: START !!!#
  for (go_timepoint_comparison in go_timepoint_comparisons) {
    
    # [2c Genes of Interest]: Loops Through the FILTERED ListComparisons From Just Above
    cat(magenta(paste('Comparison:', go_timepoint, go_timepoint_comparison,'\n'))) # Just printing to console
    geneOfInterest <- as.data.frame(results(dds, name=go_timepoint_comparison)) %>% tibble::rownames_to_column('gene_id')
    geneOfInterest <- geneOfInterest[(geneOfInterest$log2FoldChange > 2 | geneOfInterest$log2FoldChange < -2) & geneOfInterest$padj <= 0.05 & !is.na(geneOfInterest$padj), ] 
    geneOfInterest <- as.character(geneOfInterest[,'gene_id'])
    
    # [PRECAUTION]: remove any candidate genes that aren't in the go_ids dataframe from earlier (very few, <10 per 1000 on testing)
    print(paste('Removing genes not in go_ids dataframe - Before:', length(geneOfInterest))) # Just printing to console
    keep <- geneOfInterest %in% go_ids[, c('mgi_symbol')]
    geneOfInterest <- geneOfInterest[keep]
    print(paste('Removing genes not in go_ids dataframe - After:', length(geneOfInterest))) # Just printing to console
    rm(keep)
    
    # This tells TopGO at what position genes of interest appear in the 'geneUniverse' list 
    # It asks if gene of interest is in the universe and reports the answer as an integer of 0 or 1, that is then saved as a factor 0 or 1 to geneList
    # It then maps the gene ids from the Universe list as rownames so each gene id read as 0 if it is not a gene of interest, and 1 if it is
    geneList <- factor(as.integer(geneUniverse %in% geneOfInterest))
    names(geneList) <- geneUniverse

    ###################
    # Step 3: Make topGO data object
    ###################
    
    # 1) Object is topGOdata type with ontology of interest normally biological process (BP)
    # 2) allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The genes listed in this object define the gene universe.
    # 3) nodeSize: (optional) This parameter is used to prune the GO hierarchy from the terms which have < nodeSize # of annotated genes.
    # 4) annotationFun: function which maps genes identifiers to GO terms. There are a couple of annotation function included in the package trying to address the user’s needs. The annotation functions take three arguments. One of those arguments is specifying where the mappings can be found, and needs to be provided by the user. annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping.
    TopGOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = go_ids_wide)
    rm(geneList)
    
    ###################
    # Step 4: Test for significance, correct for multiple tests
    ###################
    
    # TopGo takes the GO hierarchy into account when calculating enrichment, avoiding 'inheritance' problems where general GO terms and their subterms are all significant
    # Two methods described, 'elim' and 'weight'
    # Elim discards genes with significantly enriched descendant terms
    # Weight compares significance scores of connected nodes
    # Classic tests each GO term independently
    # However: 'weight01' is the default method and is a mixture of the elim and weight methods
    TopGOdata_weight_fisher_result <- runTest(TopGOdata, algorithm='weight01', statistic='fisher') 
    
    # generate a table of results: GenTable function creates a readable summary table with the results from tests applied to the topGOdata object
    allGO <- usedGO(TopGOdata)
    all_res <- GenTable(TopGOdata, weightFisher=TopGOdata_weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))
    write.table(all_res, paste0('TopGO/noFDRadj_', go_timepoint, '_', go_timepoint_comparison, '.csv'),sep=",",quote=FALSE,row.names=FALSE)
    rm(allGO)
    
    # Generally, the p-values returned by enrichment methods in topGO are interpreted as corrected or not affected by multiple testing.
    # However, one can perform an adjustment of the p-values if they consider that it is important for their analysis.
    # Note: It is important to mention that for the methods that account for the GO topology, the problem of multiple testing is complicated.
    # Here, one computes the p-value of a GO term conditioned on the neighbouring terms. The tests are therefore not independent and the multiple testing theory does not directly apply.
    
    # fdr (Benjamini & Hochberg / BH) adjusted p-values
    # cbind p.adj to results table above
    # select only if <=0.05
    p.adj <- p.adjust(all_res$weightFisher, method="fdr")
    all_res_final <- cbind(all_res,p.adj)
    rm(all_res, p.adj)
    all_res_final <- all_res_final[order(all_res_final$p.adj),]
    all_res_final <- all_res_final[which(all_res_final$p.adj<=0.05),]
    
    # [EXPORT and reassign name for plotting]: adjp-values as .csv in TopGO folder
    write.table(all_res_final, paste0('TopGO/', go_timepoint, '_', go_timepoint_comparison, '.csv'),sep=",",quote=FALSE,row.names=FALSE)
    assign(paste0('go_result_', go_timepoint_comparison), all_res_final) # for normalised counts - therefore - any dds works
    rm(geneOfInterest, all_res_final,TopGOdata, TopGOdata_weight_fisher_result)
  }
  # End of results comparison for loop
}
# End of timepoints for loop 
rm(dds, geneUniverse, go_ids, go_ids_wide, go_timepoint, go_timepoint_comparison, go_timepoint_comparisons)


###################
# Step 5: Merge Results Tables and Plot
###################

# [5a] Join using a for loop
# This creates a list of the output dataframes from TopGO
# It uses the first one on the list to create a blank dataframe with a GO.ID column to merge with
# The for loop takes GOID, Term and the padj values and merges them into a single data frame

go_results_list <- ls(pattern = 'go_result_treat')
go_results_list
go_results <- data.frame(GO.ID = get(go_results_list[1])[,"GO.ID"]) # creates something to merge with from 1st entry

for (go_result in go_results_list) {
  temp <- get(go_result)[,c(1,2,7)]
  colnames(temp)[3] <- paste0(go_result)
  go_results <- full_join(go_results, temp)
  rm(temp)
}
rm(go_results_list, go_result)

# [5b]: Tidy go_results data.frame for complex heatmap
go_results_tidy <- as.matrix(go_results) # convert to matrix
rownames(go_results_tidy) <- paste(go_results_tidy[,'GO.ID'], go_results_tidy[,'Term']) # Merge GO.ID and Term column to use as matrix rowname
go_results_tidy <- go_results_tidy[,c(-1,-2)] # Remove GO.ID and Term columns
str(go_results_tidy) # This is a matrix full of character strings not numbers - Thanks to however TopGo saves it's results!
class(go_results_tidy) <- 'numeric' # Convert to numeric

# [5c]: MANUAL - Need to tidy up column names
colnames(go_results_tidy) <- gsub("go_result_treatment_", "", colnames(go_results_tidy)) # cleanup the length
colnames(go_results_tidy)
colnames(go_results_tidy) <- c('aPD1 vs control contralateral',
                               'aPD1 vs control injected',
                               'RP1 vs control contralateral',
                               'RP1 vs control injected',
                               'RP1aPD1 vs control contralateral',
                               'RP1aPD1 vs control injected'
                               )
# Reorder the columns (if necessary) for illustrative purposes
go_results_tidy2 <- go_results_tidy[ , c('aPD1 vs control injected',
                                         'RP1 vs control injected',
                                         'RP1aPD1 vs control injected',
                                         'aPD1 vs control contralateral',
                                         'RP1 vs control contralateral',
                                         'RP1aPD1 vs control contralateral'
                                         )]

# Order to bunch at top from left to right
go_results_tidy2 <- go_results_tidy2[order(-go_results_tidy2[,1], 
                                           -go_results_tidy2[,2],
                                           -go_results_tidy2[,3],
                                           -go_results_tidy2[,4],
                                           -go_results_tidy2[,5],
                                           -go_results_tidy2[,6]
                                           ), ]


h <- Heatmap(go_results_tidy2,
        cluster_rows = FALSE, # MUST BE TURNED OFF
        cluster_columns = FALSE, # MUST BE TURNED OFF
        column_names_side = "top",
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        col = colorRampPalette((brewer.pal(n = 7, name = "RdYlBu")))(100), # bluered(100). #rev(brewer inverts
        width = unit(ncol(go_results_tidy)*0.3, "cm"),
        height = unit(nrow(go_results_tidy)*0.225, "cm"),
        heatmap_legend_param = list(
          title = 'adjp-value',
          title_gp = gpar(fontsize = 8, font = 2),
          labels_gp = gpar(fontsize = 8),
          grid_width = unit(.4, "cm")
        ))

pdf(paste0('TopGO/', 'TopGO_enrichment', '.pdf'), width = 4.5, height = 5.5)
draw(h, heatmap_legend_side = "left")
dev.off()

rm(list = ls(pattern = 'go_|h'))
