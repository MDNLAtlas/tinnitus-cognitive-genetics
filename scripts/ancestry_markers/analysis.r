!#/usr/bin/env Rscript
library(dplyr)

PRJ_DIR <- "/home/mdnl/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/ancestry_markers"

##########################################################################################################
### Common variants as ancestry markers in LoF GBA genes
##########################################################################################################
### Preliminary tests for common variants in GBA LoF genes
f <- paste(PRJ_DIR, "GBA_LoF_genes.txt", sep="/")
genes_dt <- read.table(f, header=T)
genes <- genes_dt$SYMBOL

### Processing NFE
files <- sapply(genes, function(g) {
  pat <- genes_dt$Gene[genes_dt$SYMBOL == g]
  if(length(pat) == 0) return(NA)
  f <- list.files(PRJ_DIR, pattern = pat, full.names = TRUE)
  if(length(f) == 0) return(NA)
  return(f)
})

files <- unlist(files)

pops <- c("African.African.American", "Admixed.American", "Ashkenazi.Jewish", "East.Asian", "Middle.Eastern", "European..Finnish.", "European..non.Finnish.", "Amish", "Remaining", "South.Asian")
pops <- c(
  "African.African.American" = "afr",
  "Admixed.American" = "amr",
  "Ashkenazi.Jewish" = "asj",
  "East.Asian" = "eas",
  "Middle.Eastern" = "mid",
  "European..Finnish." = "fin",
  "European..non.Finnish." = "nfe",
  "Amish" = "ami",
  "Remaining" = "oth",
  "South.Asian" = "sas"
)
other_pops <- pops[setdiff(names(pops), "European..non.Finnish.")]

# Calculating max AF delta for each variant and filtering those with >20% difference between NFE and any other subpopulation
nfe_dt <- as.data.frame(matrix(ncol=13, nrow=0))
for (i in seq_along(files)) {

    df <- read.csv(files[i], header=T, na.strings=c(".", "", "NA"))
    df$Gene <- genes[i]

    for (p in names(pops)) {
        df[[paste0("AF.", pops[[p]])]] <- df[[paste0("Allele.Count.", p)]] / df[[paste0("Allele.Number.", p)]]
    }
    df <- df[!is.na(df$AF.nfe) & df$AF.nfe >= 0.05, ]

    df$AF_delta_max <- apply(
        abs(df[, paste0("AF.", other_pops)] - df$AF.nfe), 1, max, na.rm = TRUE
    )

    df <- df[df$AF_delta_max > 0.1, c("Gene", "gnomAD.ID", "VEP.Annotation", paste0("AF.", pops), "AF_delta_max")]
    nfe_dt <- rbind(nfe_dt, df)
}

### Processing MCI
f <- "cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.LoF_GBA_genes.tsv"
MCI_dt <- read.table(paste(PRJ_DIR, f, sep="/"), header=T, sep="\t", na.strings=c(".", "", "NA"))
MCI_dt$gnomAD.ID <- gsub("^chr", "", paste(MCI_dt$CHROM, MCI_dt$POS, MCI_dt$REF, MCI_dt$ALT, sep="-"))
colnames(MCI_dt)[colnames(MCI_dt) == "AF"] <- "AF.MCI"

f <- "cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.LoF_GBA_genes.tsv"
nonMCI_dt <- read.table(paste(PRJ_DIR, f, sep="/"), header=T, sep="\t", na.strings=c(".", "", "NA"))
nonMCI_dt$gnomAD.ID <- gsub("^chr", "", paste(nonMCI_dt$CHROM, nonMCI_dt$POS, nonMCI_dt$REF, nonMCI_dt$ALT, sep="-"))
colnames(nonMCI_dt)[colnames(nonMCI_dt) == "AF"] <- "AF.nonMCI"

# Merging with the list of ancestry markers to remove them from the analysis
m1 <- merge(nfe_dt, MCI_dt[, c("gnomAD.ID", "AF.MCI")], by="gnomAD.ID", all.x=TRUE)
m <- merge(m1, nonMCI_dt[, c("gnomAD.ID", "AF.nonMCI")], by="gnomAD.ID", all.x=TRUE)
cm <- m[complete.cases(m[, c("AF.MCI", "AF.nonMCI")]), ]

### Build result tables
cm <- cm[, c("Gene", "gnomAD.ID", "VEP.Annotation", "AF.MCI", "AF.nonMCI", "AF.nfe", "AF_delta_max", paste0("AF.", other_pops))]
colnames(cm)[colnames(cm) == "AF_delta_max"] <- "Max AF delta (NFE vs. other pops)"
wilcox.test(cm$AF.MCI, cm$AF.nonMCI)
wilcox.test(cm$AF.MCI, cm$AF.nfe)
wilcox.test(cm$AF.nonMCI, cm$AF.nfe)
for (c in 4:ncol(cm)) {
    cm[[c]] <- round(cm[[c]], 3)
}

write.csv(cm, paste(PRJ_DIR, "MCI_nonMCI_NFE.ancestry_markers.csv", sep="/"), row.names=F)

### Extract variants in chr2:178444501-178447508
m_subset <- cm[cm$Gene=="PRKRA",]
wilcox.test(m_subset$AF.MCI, m_subset$AF.nonMCI)
wilcox.test(m_subset$AF.MCI, m_subset$AF.nfe)
wilcox.test(m_subset$AF.nonMCI, m_subset$AF.nfe)

##########################################################################################################
### Common variants as ancestry markers in missense GBA genes
##########################################################################################################
f <- paste(PRJ_DIR, "GBA_missense_genes.txt", sep="/")
genes_dt <- read.table(f, header=T)
genes <- genes_dt$SYMBOL

### Processing NFE
files <- sapply(genes, function(g) {
  pat <- genes_dt$Gene[genes_dt$SYMBOL == g]
  if(length(pat) == 0) return(NA)
  f <- list.files(PRJ_DIR, pattern = pat, full.names = TRUE)
  if(length(f) == 0) return(NA)
  return(f)
})

files <- unlist(files)

pops <- c("African.African.American", "Admixed.American", "Ashkenazi.Jewish", "East.Asian", "Middle.Eastern", "European..Finnish.", "European..non.Finnish.", "Amish", "Remaining", "South.Asian")
pops <- c(
  "African.African.American" = "afr",
  "Admixed.American" = "amr",
  "Ashkenazi.Jewish" = "asj",
  "East.Asian" = "eas",
  "Middle.Eastern" = "mid",
  "European..Finnish." = "fin",
  "European..non.Finnish." = "nfe",
  "Amish" = "ami",
  "Remaining" = "oth",
  "South.Asian" = "sas"
)
other_pops <- pops[setdiff(names(pops), "European..non.Finnish.")]

# Calculating max AF delta for each variant and filtering those with >20% difference between NFE and any other subpopulation
nfe_dt <- as.data.frame(matrix(ncol=13, nrow=0))
for (i in seq_along(files)) {

    df <- read.csv(files[i], header=T, na.strings=c(".", "", "NA"))
    df$Gene <- genes[i]

    for (p in names(pops)) {
        df[[paste0("AF.", pops[[p]])]] <- df[[paste0("Allele.Count.", p)]] / df[[paste0("Allele.Number.", p)]]
    }
    df <- df[!is.na(df$AF.nfe) & df$AF.nfe >= 0.05, ]

    df$AF_delta_max <- apply(
        abs(df[, paste0("AF.", other_pops)] - df$AF.nfe), 1, max, na.rm = TRUE
    )

    df <- df[df$AF_delta_max > 0.1, c("Gene", "gnomAD.ID", "VEP.Annotation", paste0("AF.", pops), "AF_delta_max")]
    nfe_dt <- rbind(nfe_dt, df)
}

### Processing MCI
f <- "cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.missense_GBA_genes.tsv"
MCI_dt <- read.table(paste(PRJ_DIR, f, sep="/"), header=T, sep="\t", na.strings=c(".", "", "NA"))
MCI_dt$gnomAD.ID <- gsub("^chr", "", paste(MCI_dt$CHROM, MCI_dt$POS, MCI_dt$REF, MCI_dt$ALT, sep="-"))
colnames(MCI_dt)[colnames(MCI_dt) == "AF"] <- "AF.MCI"

f <- "cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.missense_GBA_genes.tsv"
nonMCI_dt <- read.table(paste(PRJ_DIR, f, sep="/"), header=T, sep="\t", na.strings=c(".", "", "NA"))
nonMCI_dt$gnomAD.ID <- gsub("^chr", "", paste(nonMCI_dt$CHROM, nonMCI_dt$POS, nonMCI_dt$REF, nonMCI_dt$ALT, sep="-"))
colnames(nonMCI_dt)[colnames(nonMCI_dt) == "AF"] <- "AF.nonMCI"

# Merging with the list of ancestry markers to remove them from the analysis
m1 <- merge(nfe_dt, MCI_dt[, c("gnomAD.ID", "AF.MCI")], by="gnomAD.ID", all.x=TRUE)
m <- merge(m1, nonMCI_dt[, c("gnomAD.ID", "AF.nonMCI")], by="gnomAD.ID", all.x=TRUE)
cm <- m[complete.cases(m[, c("AF.MCI", "AF.nonMCI")]), ]

### Build result tables
cm <- cm[, c("Gene", "gnomAD.ID", "VEP.Annotation", "AF.MCI", "AF.nonMCI", "AF.nfe", "AF_delta_max", paste0("AF.", other_pops))]
colnames(cm)[colnames(cm) == "AF_delta_max"] <- "Max AF delta (NFE vs. other pops)"
wilcox.test(cm$AF.MCI, cm$AF.nonMCI)
wilcox.test(cm$AF.MCI, cm$AF.nfe)
wilcox.test(cm$AF.nonMCI, cm$AF.nfe)
for (c in 4:ncol(cm)) {
    cm[[c]] <- round(cm[[c]], 3)
}

write.csv(cm, paste(PRJ_DIR, "MCI_nonMCI_NFE.missense_GBA_genes.ancestry_markers.csv", sep="/"), row.names=F)
