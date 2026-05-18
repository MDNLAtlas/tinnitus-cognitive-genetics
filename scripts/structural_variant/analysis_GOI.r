# revised 2025-07-29


library(synaptome.db)
setwd("~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/")

sample_dt <- read.table("script/sample_status.tsv", header = TRUE)
sample <- sample_dt[sample_dt$CNVkit=="Yes" & sample_dt$Manta=="Yes" & sample_dt$Tiddit=="Yes", "SampleID"]

# Load gene of interest
LoF_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.LoF.GBA.v3.csv", header=T)
missense_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.missense.GBA.v3.csv", header=T)
LoF_genes <- LoF_dt$SYMBOL
missense_genes <- missense_dt$SYMBOL

# Extract synaptic genes
LoF_synap_dt <- as.data.frame(getGeneInfoByName(LoF_genes))
missense_synap_dt <- as.data.frame(getGeneInfoByName(missense_genes))

GOI <- c(unique(LoF_synap_dt[,6]), unique(missense_synap_dt[,6]))
GOI_dt <- data.frame(GBA_gene=GOI)
write.table(GOI_dt, "GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.GBA_genes.v3.txt", row.names=FALSE, quote=FALSE)

# o_2
o_2 <- data.frame(matrix(ncol=(1 + length(sample)), nrow=length(GOI)))
colnames(o_2) <- c("Gene", sample)
o_2$Gene <- GOI

for (i in sample) {
    print(i)
    dt <- read.delim(paste0("merged/", i, "/", i, ".merged.flt.annotSV.vcf.gz.tsv"), header=T, na.strings = c("NA", ""))

    # o_2
    idx <- 1
    for (g in GOI) {
        print(g)
        dt_1 <- dt[dt$Annotation_mode=="split" & dt$Gene_name == g, ]
        dt_1 <- dt_1[rowSums(is.na(dt_1)) != ncol(dt_1), ]
        if (nrow(dt_1) > 0) {
            SV_id <- unique(paste0(dt_1$AnnotSV_ID, "-", dt_1$ID))
            o_2[idx, i] <- paste(SV_id, collapse = "; ")
        }
        idx <- idx + 1
    }
}

dir.create("SV_CNV_output", showWarnings = FALSE)
write.csv(o_2, "SV_CNV_output/SV_GOI_overlap_revised.v3.csv", row.names = FALSE)

# 001 TIDDIT

# o_2
o_2 <- data.frame(matrix(ncol=(1 + length(sample)), nrow=length(GOI)))
colnames(o_2) <- c("Gene", sample)
o_2$Gene <- GOI

for (i in sample) {
    print(i)
    dt <- read.delim(paste0("merged/", i, "/", i, ".merged.flt.001.annotSV.vcf.gz.tsv"), header=T, na.strings = c("NA", ""))


    # o_2
    idx <- 1
    for (g in GOI) {
        print(g)
        dt_1 <- dt[dt$Annotation_mode=="split" & dt$Gene_name == g, ]
        dt_1 <- dt_1[rowSums(is.na(dt_1)) != ncol(dt_1), ]
        if (nrow(dt_1) > 0) {
            SV_id <- unique(paste0(dt_1$AnnotSV_ID, "-", dt_1$ID))
            o_2[idx, i] <- paste(SV_id, collapse = "; ")
        }
        idx <- idx + 1
    }
}

row_vec <- as.character(o_2[1, ])
split_items <- unlist(strsplit(row_vec, "; "))
result <- sapply(split_items, function(s) {
  parts <- strsplit(s, "_")[[1]][1:4]
  paste(parts, collapse = "_")
})
table(result)

dir.create("SV_CNV_output", showWarnings = FALSE)
write.csv(o_2, "SV_CNV_output/SV_GOI_overlap_revised.TIDDIT.001.csv", row.names = FALSE)

# 010 Manta

# o_2
o_2 <- data.frame(matrix(ncol=(1 + length(sample)), nrow=length(GOI)))
colnames(o_2) <- c("Gene", sample)
o_2$Gene <- GOI

for (i in sample) {
    print(i)
    dt <- read.delim(paste0("merged/", i, "/", i, ".merged.flt.010.annotSV.vcf.gz.tsv"), header=T, na.strings = c("NA", ""))

    # o_2
    idx <- 1
    for (g in GOI) {
        print(g)
        dt_1 <- dt[dt$Annotation_mode=="split" & dt$Gene_name == g, ]
        dt_1 <- dt_1[rowSums(is.na(dt_1)) != ncol(dt_1), ]
        if (nrow(dt_1) > 0) {
            SV_id <- unique(paste0(dt_1$AnnotSV_ID, "-", dt_1$ID))
            o_2[idx, i] <- paste(SV_id, collapse = "; ")
        }
        idx <- idx + 1
    }
}

row_vec <- as.character(o_2[10, ])
split_items <- unlist(strsplit(row_vec, "; "))
result <- sapply(split_items, function(s) {
  parts <- strsplit(s, "_")[[1]][1:4]
  paste(parts, collapse = "_")
})
table(result)

dir.create("SV_CNV_output", showWarnings = FALSE)
write.csv(o_2, "SV_CNV_output/SV_GOI_overlap_revised.Manta.010.csv", row.names = FALSE)

# cnvkit
o_2 <- data.frame(matrix(ncol=(1 + length(sample)), nrow=length(GOI)))
colnames(o_2) <- c("Gene", sample)
o_2$Gene <- GOI

for (i in sample) {
    print(i)
    dt <- read.delim(paste0("merged/", i, "/", i, ".cnvcall.withID.annotSV.vcf.gz.tsv"), header=T, na.strings = c("NA", ""))

    # o_2
    idx <- 1
    for (g in GOI) {
        print(g)
        dt_1 <- dt[dt$Annotation_mode=="split" & dt$Gene_name == g, ]
        dt_1 <- dt_1[rowSums(is.na(dt_1)) != ncol(dt_1), ]
        if (nrow(dt_1) > 0) {
            SV_id <- unique(paste0(dt_1$AnnotSV_ID, "-", dt_1$ID))
            o_2[idx, i] <- paste(SV_id, collapse = "; ")
        }
        idx <- idx + 1
    }
}

dir.create("SV_CNV_output", showWarnings = FALSE)
write.csv(o_2, "SV_CNV_output/CNV_GOI_overlap.csv", row.names = FALSE)
