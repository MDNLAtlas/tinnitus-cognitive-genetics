.libPaths("/home/mdnl/R/x86_64-pc-linux-gnu-library/4.4")

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(VennDiagram)
library(tidyr)
library(dplyr)
library(gridExtra)

setwd("/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1")

### Define a function to convert entrezID to UCSC gene symbol
convert_entrez_to_symbol <- function(entrez_ids) {
    entrez_list <- unlist(strsplit(entrez_ids, "/"))
    
    # Use bitr to convert Entrez IDs to Gene Symbols
    gene_symbols <- bitr(entrez_list, 
                         fromType = "ENTREZID", 
                         toType = "SYMBOL", 
                         OrgDb = org.Hs.eg.db)
    
    return(paste(gene_symbols$SYMBOL, collapse = "/"))
}

### Convert geneRatio and bgRatio into numeric
ratioStr2ratioNum <- function(ratioStr) {
    num <- as.numeric(unlist(strsplit(ratioStr, "/")))
    return(num[1] / num[2])
}

LoF_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.LoF.GBA.v3.csv", header=T)
missense_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.missense.GBA.v3.csv", header=T)

genes <- unique(c(LoF_dt$SYMBOL, missense_dt$SYMBOL))
entrez_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

### GO BP pathway
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2,
                            pAdjustMethod = "BH")

go_results_df <- as.data.frame(go_results)
go_results_df$geneSymbol <- sapply(go_results_df$geneID,
    function(entrezID_list) convert_entrez_to_symbol(entrezID_list))
geneRatio <- sapply(go_results_df$GeneRatio, ratioStr2ratioNum)
BgRatio <- sapply(go_results_df$BgRatio, ratioStr2ratioNum)
go_results_df$OR <- (geneRatio / (1-geneRatio)) / (BgRatio / (1-BgRatio))
write.csv(go_results_df, "GO_MF.csv")

library(ggplot2)
library(dplyr)
library(forcats)

# Prepare GO data (example: top 20 terms)
go_plot_df <- go_results_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(Description = fct_reorder(Description, GeneRatio),
    GeneRatio = sapply(GeneRatio, function(x) {
        parts <- as.numeric(unlist(strsplit(x, "/")))
        parts[1] / parts[2]
    }),
    
    # Ensure Count is integer
    Count = as.integer(Count))  # reorder by GeneRatio

# Draw bubble plot
png("figures/GO_MF.png", width=10*300, height=5*300, res=300)
ggplot(go_plot_df, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. P value") +
  scale_size(range = c(3, 10)) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    title = "GO Molecular Function (MF) Enrichment",
    size = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()



## KEGG pathway
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                                organism = "hsa",  # for human
                                pvalueCutoff = 0.05)
kegg_results_df <- as.data.frame(kegg_results)
kegg_results_df$geneSymbol <- sapply(kegg_results_df$geneID,
    function(entrezID_list) convert_entrez_to_symbol(entrezID_list)) 
geneRatio <- sapply(kegg_results_df$GeneRatio, ratioStr2ratioNum)
BgRatio <- sapply(kegg_results_df$BgRatio, ratioStr2ratioNum)
kegg_results_df$OR <- (geneRatio / (1-geneRatio)) / (BgRatio / (1-BgRatio))