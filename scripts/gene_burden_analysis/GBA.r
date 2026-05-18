# revised 20/08/2025
library(tidyr)
library(dplyr)

# SNV
setwd("/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1")

cons_type <- c("missense", "LoF")

for (c in cons_type) {
    print(c)
    df <- read.table(paste0("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
    df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
    # check if any duplicates
    # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
    
    ## GBA
    # df.gba <- df
    # length(unique(df$SYMBOL))
    # 417 for LoF

    df.gba <- df %>%
        mutate(n_all = n_het + n_hom) %>%
        group_by(Gene, SYMBOL) %>%
        filter(!all(n_all == 1)) %>%
        ungroup()

    df.gba <- df.gba %>%
        group_by(Gene, SYMBOL) %>%
        summarise(No.of.variants=n(),
                    Pos=paste(Pos, collapse=";"),
                    Consequence=paste(Consequence, collapse=";"),
                    EXON=paste(EXON, collapse=";"),
                    HGVSc=paste(HGVSc, collapse=";"),
                    HGVSp=paste(HGVSp, collapse=";"),
                    CDS_position=paste(CDS_position, collapse=";"),
                    Amino_acids=paste(Amino_acids, collapse=";"),
                    Codons=paste(Codons, collapse=";"),
                    cDNA_position=paste(cDNA_position, collapse=";"),
                    Protein_position=paste(Protein_position, collapse=";"),
                    CANONICAL=paste(CANONICAL, collapse=";"),
                    HGNC_ID=paste(HGNC_ID, collapse=";"),
                    pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                    CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                    gnomADg=paste(gnomADg, collapse=";"),
                    n_het=paste(n_het, collapse=";"),
                    n_hom=paste(n_hom, collapse=";"),
                    samples_het=paste(samples_het, collapse=";"),
                    samples_hom=paste(samples_hom, collapse=";"),
                    gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                    gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                    AC=sum(AC),
                    AN=sum(AN),
                    RC=sum(AN) - sum(AC),
                    gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                    gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                    gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                    gnomADg_AC=sum(gnomADg_AC),
                    gnomADg_AN=sum(gnomADg_AN),
                    gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
    df.gba <- as.data.frame(df.gba)
    # m <- match(df.gba$Gene, df$Gene)
    # df.gba$SYMBOL <- df$SYMBOL[m]

    ### NFE gnomeADg
    t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
        t_res <- matrix(x, nrow = 2, byrow=T)
        OR <- (x[1] * x[4]) / (x[2] * x[3])
        logOR <- log(OR)
        SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
        z <- logOR / SE_logOR
        p_value <- 2 * pnorm(-abs(z))
        CI_lower <- exp(logOR - 1.96 * SE_logOR)
        CI_upper <- exp(logOR + 1.96 * SE_logOR)
        return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
    })
    t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
    adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
    EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

    ### ALL gnomeADg
    t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
        t_res <- matrix(x, nrow = 2, byrow=T)
        OR <- (x[1] * x[4]) / (x[2] * x[3])
        logOR <- log(OR)
        SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
        z <- logOR / SE_logOR
        p_value <- 2 * pnorm(-abs(z))
        CI_lower <- exp(logOR - 1.96 * SE_logOR)
        CI_upper <- exp(logOR + 1.96 * SE_logOR)
        return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
    })
    t.gnomADg <- as.data.frame(t(t.gnomADg))
    adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
    EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

    df.res <- cbind(df.gba, 
                    p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                    FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                    OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                    CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                    CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                    EF_NFE_gnomADg=EF.gnomADg.nfe,   
                    p_value_ALL_gnomADg=t.gnomADg$pVal, 
                    FDR_ALL_gnomADg=adj.pVal.gnomADg,
                    OR_ALL_gnomADg=t.gnomADg$OR,
                    CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                    CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                    EF_ALL_gnomADg=EF.gnomADg,
                    IS.SIG=rep("ns", nrow(df.gba)))

    df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
    df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
    df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

    df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
    write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.", c, ".GBA.v3.csv"))

}

# GBA_GUEF
for (q in c("Q1", "Q4")) {
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.csv"))

    }

}

# GBA_THI
for (q in c(Q1, Q4)) {
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.csv"))

    }
}

# GBA THI (5 stages)
for (q in paste0("g", 1:5)) {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        count_unique_samples <- function(samples_het, samples_hom) {
            combined <- paste(samples_het, samples_hom, sep = ";")
            samples <- unlist(strsplit(combined, ";|\\|"))
            samples <- samples[!samples %in% c("NA", "", " ")]
            return(length(unique(samples)))
        }
        df.res$sample_count <- mapply(count_unique_samples, df.res$samples_het, df.res$samples_hom)
  
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        df.out <- cbind(df.res[, c("SYMBOL", "No.of.variants", "sample_count", "IS.SIG")],
                        `OR_gnomADg_NFE`=paste0(round(df.res$OR_NFE_gnomADg,2), "\n(", round(df.res$CI_lower_NFE_gnomADg,2), " - ", round(df.res$CI_upper_NFE_gnomADg,2), ")"),
                        `FDR_gnomADg_NFE`=sprintf("%.2E", df.res$FDR_NFE_gnomADg),
                        `EF_gnomADg_NFE`=round(df.res$EF_NFE_gnomADg,2),
                        `OR_gnomADg_ALL`=paste0(round(df.res$OR_ALL_gnomADg,2), "\n(", round(df.res$CI_lower_ALL_gnomADg,2), " - ", round(df.res$CI_upper_ALL_gnomADg,2), ")"),
                        `FDR_gnomADg_ALL`=sprintf("%.2E", df.res$FDR_ALL_gnomADg), 
                        `EF_gnomADg_ALL`=round(df.res$EF_ALL_gnomADg,2))
        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.csv"))
        write.csv(df.out[df.out$IS.SIG %in% c("NFE", "NFE/ALL"), -4], paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.formated.csv"), row.names = F)

    }

}

# GBA THI (5 stages)
for (q in c("mild", "severe")) {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        count_unique_samples <- function(samples_het, samples_hom) {
            combined <- paste(samples_het, samples_hom, sep = ";")
            samples <- unlist(strsplit(combined, ";|\\|"))
            samples <- samples[!samples %in% c("NA", "", " ")]
            return(length(unique(samples)))
        }
        df.res$sample_count <- mapply(count_unique_samples, df.res$samples_het, df.res$samples_hom)
  
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        df.out <- cbind(df.res[, c("SYMBOL", "No.of.variants", "sample_count", "IS.SIG")],
                        `OR_gnomADg_NFE`=paste0(round(df.res$OR_NFE_gnomADg,2), "\n(", round(df.res$CI_lower_NFE_gnomADg,2), " - ", round(df.res$CI_upper_NFE_gnomADg,2), ")"),
                        `FDR_gnomADg_NFE`=sprintf("%.2E", df.res$FDR_NFE_gnomADg),
                        `EF_gnomADg_NFE`=round(df.res$EF_NFE_gnomADg,2),
                        `OR_gnomADg_ALL`=paste0(round(df.res$OR_ALL_gnomADg,2), "\n(", round(df.res$CI_lower_ALL_gnomADg,2), " - ", round(df.res$CI_upper_ALL_gnomADg,2), ")"),
                        `FDR_gnomADg_ALL`=sprintf("%.2E", df.res$FDR_ALL_gnomADg), 
                        `EF_gnomADg_ALL`=round(df.res$EF_ALL_gnomADg,2))
        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.csv"))
        write.csv(df.out[df.out$IS.SIG %in% c("NFE", "NFE/ALL"), -4], paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.formated.csv"), row.names = F)

    }

}

# GBA_GUEF 4 stages
for (q in paste0("g", 1:4)) {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        count_unique_samples <- function(samples_het, samples_hom) {
            combined <- paste(samples_het, samples_hom, sep = ";")
            samples <- unlist(strsplit(combined, ";|\\|"))
            samples <- samples[!samples %in% c("NA", "", " ")]
            return(length(unique(samples)))
        }
        df.res$sample_count <- mapply(count_unique_samples, df.res$samples_het, df.res$samples_hom)
  
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        df.out <- cbind(df.res[, c("SYMBOL", "No.of.variants", "sample_count", "IS.SIG")],
                        `OR_gnomADg_NFE`=paste0(round(df.res$OR_NFE_gnomADg,2), "\n(", round(df.res$CI_lower_NFE_gnomADg,2), " - ", round(df.res$CI_upper_NFE_gnomADg,2), ")"),
                        `FDR_gnomADg_NFE`=sprintf("%.2E", df.res$FDR_NFE_gnomADg),
                        `EF_gnomADg_NFE`=round(df.res$EF_NFE_gnomADg,2),
                        `OR_gnomADg_ALL`=paste0(round(df.res$OR_ALL_gnomADg,2), "\n(", round(df.res$CI_lower_ALL_gnomADg,2), " - ", round(df.res$CI_upper_ALL_gnomADg,2), ")"),
                        `FDR_gnomADg_ALL`=sprintf("%.2E", df.res$FDR_ALL_gnomADg), 
                        `EF_gnomADg_ALL`=round(df.res$EF_ALL_gnomADg,2))

        # write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.csv"))
        write.csv(df.out[df.out$IS.SIG %in% c("NFE", "NFE/ALL"), -4], paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.formated.csv"), row.names = F)

    }

}

for (q in c("severe")) {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        count_unique_samples <- function(samples_het, samples_hom) {
            combined <- paste(samples_het, samples_hom, sep = ";")
            samples <- unlist(strsplit(combined, ";|\\|"))
            samples <- samples[!samples %in% c("NA", "", " ")]
            return(length(unique(samples)))
        }
        df.res$sample_count <- mapply(count_unique_samples, df.res$samples_het, df.res$samples_hom)
  
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        df.out <- cbind(df.res[, c("SYMBOL", "No.of.variants", "sample_count", "IS.SIG")],
                        `OR_gnomADg_NFE`=paste0(round(df.res$OR_NFE_gnomADg,2), "\n(", round(df.res$CI_lower_NFE_gnomADg,2), " - ", round(df.res$CI_upper_NFE_gnomADg,2), ")"),
                        `FDR_gnomADg_NFE`=sprintf("%.2E", df.res$FDR_NFE_gnomADg),
                        `EF_gnomADg_NFE`=round(df.res$EF_NFE_gnomADg,2),
                        `OR_gnomADg_ALL`=paste0(round(df.res$OR_ALL_gnomADg,2), "\n(", round(df.res$CI_lower_ALL_gnomADg,2), " - ", round(df.res$CI_upper_ALL_gnomADg,2), ")"),
                        `FDR_gnomADg_ALL`=sprintf("%.2E", df.res$FDR_ALL_gnomADg), 
                        `EF_gnomADg_ALL`=round(df.res$EF_ALL_gnomADg,2))

        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.csv"))
        write.csv(df.out[df.out$IS.SIG %in% c("NFE", "NFE/ALL"), -4], paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.formated.csv"), row.names = F)

    }

}

# PHQ9
for (q in paste0("g", 1:4)) {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        count_unique_samples <- function(samples_het, samples_hom) {
            combined <- paste(samples_het, samples_hom, sep = ";")
            samples <- unlist(strsplit(combined, ";|\\|"))
            samples <- samples[!samples %in% c("NA", "", " ")]
            return(length(unique(samples)))
        }
        df.res$sample_count <- mapply(count_unique_samples, df.res$samples_het, df.res$samples_hom)
  
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        df.out <- cbind(df.res[, c("SYMBOL", "No.of.variants", "sample_count", "IS.SIG")],
                        `OR_gnomADg_NFE`=paste0(round(df.res$OR_NFE_gnomADg,2), "\n(", round(df.res$CI_lower_NFE_gnomADg,2), " - ", round(df.res$CI_upper_NFE_gnomADg,2), ")"),
                        `FDR_gnomADg_NFE`=sprintf("%.2E", df.res$FDR_NFE_gnomADg),
                        `EF_gnomADg_NFE`=round(df.res$EF_NFE_gnomADg,2),
                        `OR_gnomADg_ALL`=paste0(round(df.res$OR_ALL_gnomADg,2), "\n(", round(df.res$CI_lower_ALL_gnomADg,2), " - ", round(df.res$CI_upper_ALL_gnomADg,2), ")"),
                        `FDR_gnomADg_ALL`=sprintf("%.2E", df.res$FDR_ALL_gnomADg), 
                        `EF_gnomADg_ALL`=round(df.res$EF_ALL_gnomADg,2))

        # write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_", q, ".", c, ".GBA.csv"))
        write.csv(df.out[df.out$IS.SIG %in% c("NFE", "NFE/ALL"), -4], paste0("GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_", q, ".", c, ".GBA.formated.csv"), row.names = F)
    }

}

# MOCA
for (q in "normal") {
    print(q)
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df$AC_Hom <- sapply(strsplit(as.character(df$AC_Hom), ","), function(vec) {
        paste(as.numeric(vec) / 2, collapse = ",")
        })
        colnames(df)[colnames(df) == "AC_Het"] <- "n_het"
        colnames(df)[colnames(df) == "AC_Hom"] <- "n_hom"

        df$n_het <- as.numeric(df$n_het)
        df$n_hom <- as.numeric(df$n_hom)

        samples <- colnames(df)[57:ncol(df)]

        get_sample_groups <- function(row) {
            het_samples <- c()
            hom_samples <- c()
        
            for (a in 1:length(strsplit(row[["ALT"]], ",")[[1]])) {

                hom_s <- c()
                het_s <- c()
                for (s in samples) {
                    gt <- row[[s]]
                    # s <- gsub("\\.", "-", s)
                    if (is.na(gt) || gt %in% c("0/0", "0|0", "./.", ".|.")) next
                    
                    alleles <- unlist(strsplit(gt, "[/|]"))
                    if (length(alleles) != 2 || any(alleles == ".")) next
                    
                    if (alleles[1] == alleles[2] && alleles[1] == a) {
                        hom_s <- c(hom_s, s)
                    } else if ((alleles[1] == a || alleles[2] == a) && alleles[1] != alleles[2]) {
                        het_s <- c(het_s, s)
                    }
                }

                if (length(het_s) == 0) {
                    het_s <- NA
                } else {
                    het_s <- paste(het_s, collapse = "|")
                }

                if (length(hom_s) == 0) {
                    hom_s <- NA
                } else {    
                    hom_s <- paste(hom_s, collapse = "|")
                }

                het_samples <- c(het_samples, het_s)
                hom_samples <- c(hom_samples, hom_s)

            }

            het_samples <- paste(het_samples, collapse = ",")
            hom_samples <- paste(hom_samples, collapse = ",")

            return(list(samples_het = het_samples, samples_hom = hom_samples))
        }

        result <- apply(df, 1, get_sample_groups)
        df$samples_het <- sapply(result, function(x) x$samples_het)
        df$samples_hom <- sapply(result, function(x) x$samples_hom)

        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        # check if any duplicates
        # dim(df)[1] == dim(unique(df[,c(9,56)]))[1] # based on SYMBOL and Pos
        
        ## GBA
        # df.gba <- df
        # length(unique(df$SYMBOL))
        # 417 for LoF

        df.gba <- df %>%
            mutate(n_all = n_het + n_hom) %>%
            group_by(Gene, SYMBOL) %>%
            filter(!all(n_all == 1)) %>%
            ungroup()

        df.gba <- df.gba %>%
            group_by(Gene, SYMBOL) %>%
            summarise(No.of.variants=n(),
                        Pos=paste(Pos, collapse=";"),
                        Consequence=paste(Consequence, collapse=";"),
                        EXON=paste(EXON, collapse=";"),
                        HGVSc=paste(HGVSc, collapse=";"),
                        HGVSp=paste(HGVSp, collapse=";"),
                        CDS_position=paste(CDS_position, collapse=";"),
                        Amino_acids=paste(Amino_acids, collapse=";"),
                        Codons=paste(Codons, collapse=";"),
                        cDNA_position=paste(cDNA_position, collapse=";"),
                        Protein_position=paste(Protein_position, collapse=";"),
                        CANONICAL=paste(CANONICAL, collapse=";"),
                        HGNC_ID=paste(HGNC_ID, collapse=";"),
                        pLI_gene_value=paste(pLI_gene_value, collapse=";"),
                        CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                        gnomADg=paste(gnomADg, collapse=";"),
                        n_het=paste(n_het, collapse=";"),
                        n_hom=paste(n_hom, collapse=";"),
                        samples_het=paste(samples_het, collapse=";"),
                        samples_hom=paste(samples_hom, collapse=";"),
                        gnomADg_AF_nfe=paste(gnomADg_AF_nfe, collapse=";"),
                        gnomADg_AF=paste(gnomADg_AF, collapse=";"),
                        AC=sum(AC),
                        AN=sum(AN),
                        RC=sum(AN) - sum(AC),
                        gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                        gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                        gnomADg_RC_nfe=sum(gnomADg_AN_nfe) - sum(gnomADg_AC_nfe),
                        gnomADg_AC=sum(gnomADg_AC),
                        gnomADg_AN=sum(gnomADg_AN),
                        gnomADg_RC=sum(gnomADg_AN) - sum(gnomADg_AC), .groups="drop")
        df.gba <- as.data.frame(df.gba)
        # m <- match(df.gba$Gene, df$Gene)
        # df.gba$SYMBOL <- df$SYMBOL[m]

        ### NFE gnomeADg
        t.gnomADg.nfe <- apply(df.gba[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
        adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
        EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

        ### ALL gnomeADg
        t.gnomADg <- apply(df.gba[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
            t_res <- matrix(x, nrow = 2, byrow=T)
            OR <- (x[1] * x[4]) / (x[2] * x[3])
            logOR <- log(OR)
            SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            z <- logOR / SE_logOR
            p_value <- 2 * pnorm(-abs(z))
            CI_lower <- exp(logOR - 1.96 * SE_logOR)
            CI_upper <- exp(logOR + 1.96 * SE_logOR)
            return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
        })
        t.gnomADg <- as.data.frame(t(t.gnomADg))
        adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
        EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

        df.res <- cbind(df.gba, 
                        p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                        FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                        OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                        CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                        CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                        EF_NFE_gnomADg=EF.gnomADg.nfe,   
                        p_value_ALL_gnomADg=t.gnomADg$pVal, 
                        FDR_ALL_gnomADg=adj.pVal.gnomADg,
                        OR_ALL_gnomADg=t.gnomADg$OR,
                        CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                        CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                        EF_ALL_gnomADg=EF.gnomADg,
                        IS.SIG=rep("ns", nrow(df.gba)))

        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
        df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
        df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

        df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
        print(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),]$SYMBOL)
        write.csv(df.res[df.res$IS.SIG %in% c("NFE", "NFE/ALL"),], paste0("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_", q, ".", c, ".GBA.csv"))

    }

}

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(cowplot)
# setwd("/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1")

# MCI vs. non-MCI
n.MCI <- 75
n.nonMCI <- 201

sv <- list()
sv[["Q1"]] <- data.frame(
    SYMBOL = c("PRKRA"),
    Pos = c("chr2:178444501-178447508:DEL"),
    AC = c(14),
    AN = c(n.MCI * 2),
    gnomADg_AC_nfe = c(10),
    gnomADg_AN_nfe = c(49918),
    gnomADg_AC = c(37),
    gnomADg_AN = c(108658)
)

sv[["normal"]] <- data.frame(
    SYMBOL = c("PRKRA"),
    Pos = c("chr2:178444501-178447508:DEL"),
    AC = c(6),
    AN = c(n.nonMCI * 2),
    gnomADg_AC_nfe = c(10),
    gnomADg_AN_nfe = c(49918),
    gnomADg_AC = c(37),
    gnomADg_AN = c(108658)
)

synap_gba <- read.table("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.GBA_genes.v3.txt", header=F)
cons_type <- c("LoF", "missense")

df.sva <- list()
df.gba <- list()
for (q in c("Q1", "normal")) {
    print(q)
    df.full <- data.frame()
    for (c in cons_type) {
        print(c)
        df <- read.table(paste0("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_", q, ".", c, ".tsv"), header = T, sep = "\t", na.strings=c(".", "", "NA", "N/A"))
        df <- df[df$SYMBOL %in% synap_gba$V1, ]
        df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)
        df.full <- rbind(df.full, df[, c("SYMBOL", "Pos", "AC", "AN", "gnomADg_AC_nfe", "gnomADg_AN_nfe", "gnomADg_AC", "gnomADg_AN")])
    }

    df.full <- rbind(df.full, sv[[q]]) 
    df.full$RC <- df.full$AN - df.full$AC
    df.full$gnomADg_RC_nfe <- df.full$gnomADg_AN_nfe - df.full$gnomADg_AC_nfe
    df.full$gnomADg_RC <- df.full$gnomADg_AN - df.full$gnomADg_AC

    df.sva[[q]] <- df.full
    df.gba[[q]] <- df.full %>%
        group_by(SYMBOL) %>%
        summarise(No.of.variants=n(),
                    Pos=paste(Pos, collapse=";"),
                    AC=sum(AC),
                    AN=sum(AN),
                    RC=sum(RC),
                    gnomADg_AC_nfe=sum(gnomADg_AC_nfe),
                    gnomADg_AN_nfe=sum(gnomADg_AN_nfe),
                    gnomADg_RC_nfe=sum(gnomADg_RC_nfe),
                    gnomADg_AC=sum(gnomADg_AC),
                    gnomADg_AN=sum(gnomADg_AN),
                    gnomADg_RC=sum(gnomADg_RC), .groups="drop")
    df.gba[[q]] <- as.data.frame(df.gba[[q]])
}

df <- list()
df[["gba"]] <- left_join(df.gba[["Q1"]], df.gba[["normal"]], by=c("SYMBOL"), suffix=c(".MCI", ".nonMCI"))
df[["sva"]] <- left_join(df.sva[["Q1"]], df.sva[["normal"]], by=c("Pos"), suffix=c(".MCI", ".nonMCI"))

for (i in names(df)) {
    df1 <- df[[i]]

    if (i == "gba") {
        df1 <- df1[, -which(colnames(df1) %in% c("No.of.variants.nonMCI"))]
        df1$Pos.nonMCI <- ifelse(is.na(df1$Pos.nonMCI), df1$Pos.MCI, df1$Pos.nonMCI)
    } else {
        df1$SYMBOL.nonMCI <- ifelse(is.na(df1$SYMBOL.nonMCI), df1$SYMBOL.MCI, df1$SYMBOL.nonMCI)
    }
    df1$AC.nonMCI <- ifelse(is.na(df1$AC.nonMCI), 0, df1$AC.nonMCI)
    df1$AN.nonMCI <- ifelse(is.na(df1$AN.nonMCI), df1$AN.MCI, df1$AN.nonMCI)
    df1$RC.nonMCI <- df1$AN.nonMCI - df1$AC.nonMCI
    df1$gnomADg_AC_nfe.nonMCI <- ifelse(is.na(df1$gnomADg_AC_nfe.nonMCI), df1$gnomADg_AC_nfe.MCI, df1$gnomADg_AC_nfe.nonMCI)
    df1$gnomADg_AN_nfe.nonMCI <- ifelse(is.na(df1$gnomADg_AN_nfe.nonMCI), df1$gnomADg_AN_nfe.MCI, df1$gnomADg_AN_nfe.nonMCI)
    df1$gnomADg_RC_nfe.nonMCI <- df1$gnomADg_AN_nfe.nonMCI - df1$gnomADg_AC_nfe.nonMCI
    df1$gnomADg_AC.nonMCI <- ifelse(is.na(df1$gnomADg_AC.nonMCI), df1$gnomADg_AC.MCI, df1$gnomADg_AC.nonMCI)
    df1$gnomADg_AN.nonMCI <- ifelse(is.na(df1$gnomADg_AN.nonMCI), df1$gnomADg_AN.MCI, df1$gnomADg_AN.nonMCI)
    df1$gnomADg_RC.nonMCI <- df1$gnomADg_AN.nonMCI - df1$gnomADg_AC.nonMCI

    ### MCI vs. non-MCI association tests
    t <- apply(df1[, c("AC.MCI", "RC.MCI", "AC.nonMCI", "RC.nonMCI")], 1, function(x) {
        t_res <- matrix(c(x["AC.MCI"], x["RC.MCI"], x["AC.nonMCI"], x["RC.nonMCI"]), nrow = 2, byrow=T)
        OR <- (x["AC.MCI"] * x["RC.nonMCI"]) / (x["RC.MCI"] * x["AC.nonMCI"])
        logOR <- log(OR)
        SE_logOR <- sqrt(1/x["AC.MCI"] + 1/x["RC.MCI"] + 1/x["AC.nonMCI"] + 1/x["RC.nonMCI"])
        z <- logOR / SE_logOR
        p_value <- 2 * pnorm(-abs(z))
        CI_lower <- exp(logOR - 1.96 * SE_logOR)
        CI_upper <- exp(logOR + 1.96 * SE_logOR)
        return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
    })
    t <- as.data.frame(t(t))
    adj.pVal <- p.adjust(t$pVal, method = "bonferroni")
    EF <- (t$OR-1)/t$OR

    ### NFE gnomeADg MCI
    t.gnomADg.nfe.MCI <- apply(df1[, c("AC.MCI", "RC.MCI", "gnomADg_AC_nfe.MCI", "gnomADg_RC_nfe.MCI")], 1, function(x) {
        t_res <- matrix(x, nrow = 2, byrow=T)
        OR <- (x[1] * x[4]) / (x[2] * x[3])
        logOR <- log(OR)
        SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
        z <- logOR / SE_logOR
        p_value <- 2 * pnorm(-abs(z))
        CI_lower <- exp(logOR - 1.96 * SE_logOR)
        CI_upper <- exp(logOR + 1.96 * SE_logOR)
        return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper, se=SE_logOR))
    })
    t.gnomADg.nfe.MCI <- as.data.frame(t(t.gnomADg.nfe.MCI))
    adj.pVal.gnomADg.nfe.MCI <- p.adjust(t.gnomADg.nfe.MCI$pVal, method = "bonferroni")
    EF.gnomADg.nfe.MCI <- (t.gnomADg.nfe.MCI$OR-1)/t.gnomADg.nfe.MCI$OR

    ### NFE gnomeADg nonMCI
    t.gnomADg.nfe.nonMCI <- apply(df1[, c("AC.nonMCI", "RC.nonMCI", "gnomADg_AC_nfe.nonMCI", "gnomADg_RC_nfe.nonMCI")], 1, function(x) {
        t_res <- matrix(x, nrow = 2, byrow=T)
        OR <- (x[1] * x[4]) / (x[2] * x[3])
        logOR <- log(OR)
        SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
        z <- logOR / SE_logOR
        p_value <- 2 * pnorm(-abs(z))
        CI_lower <- exp(logOR - 1.96 * SE_logOR)
        CI_upper <- exp(logOR + 1.96 * SE_logOR)
        return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper, se=SE_logOR))
    })
    t.gnomADg.nfe.nonMCI <- as.data.frame(t(t.gnomADg.nfe.nonMCI))
    adj.pVal.gnomADg.nfe.nonMCI <- p.adjust(t.gnomADg.nfe.nonMCI$pVal, method = "bonferroni")
    EF.gnomADg.nfe.nonMCI <- (t.gnomADg.nfe.nonMCI$OR-1)/t.gnomADg.nfe.nonMCI$OR

    z <- (log(t.gnomADg.nfe.MCI$OR) - log(t.gnomADg.nfe.nonMCI$OR)) / sqrt(t.gnomADg.nfe.MCI$se^2 + t.gnomADg.nfe.nonMCI$se^2)
    p_diff <- 2 * pnorm(-abs(z))
    adj.p_diff <- p.adjust(p_diff, method = "BH") 


    df[[i]] <- cbind(df1, 
                        p_value_NFE_gnomADg.MCI=t.gnomADg.nfe.MCI$pVal, 
                        FDR_NFE_gnomADg.MCI=adj.pVal.gnomADg.nfe.MCI, 
                        OR_NFE_gnomADg.MCI=t.gnomADg.nfe.MCI$OR,
                        CI_lower_NFE_gnomADg.MCI=t.gnomADg.nfe.MCI$CI_lower,
                        CI_upper_NFE_gnomADg.MCI=t.gnomADg.nfe.MCI$CI_upper,
                        EF_NFE_gnomADg.MCI=EF.gnomADg.nfe.MCI,   
                        p_value_NFE_gnomADg.nonMCI=t.gnomADg.nfe.nonMCI$pVal, 
                        FDR_NFE_gnomADg.nonMCI=adj.pVal.gnomADg.nfe.nonMCI, 
                        OR_NFE_gnomADg.nonMCI=t.gnomADg.nfe.nonMCI$OR,
                        CI_lower_NFE_gnomADg.nonMCI=t.gnomADg.nfe.nonMCI$CI_lower,
                        CI_upper_NFE_gnomADg.nonMCI=t.gnomADg.nfe.nonMCI$CI_upper,
                        EF_NFE_gnomADg.nonMCI=EF.gnomADg.nfe.nonMCI,   
                        p_value=t$pVal, 
                        FDR=adj.pVal, 
                        OR=t$OR,
                        CI_lower=t$CI_lower,
                        CI_upper=t$CI_upper,
                        EF=EF,
                        p_diff=p_diff,
                        FDR_diff=adj.p_diff)
}

write.csv(df[["sva"]], "GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MCIvsnonMCI.SVA.csv", row.names=F)
write.csv(df[["gba"]], "GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MCIvsnonMCI.GBA.csv", row.names=F)

# plotting OR
p <- list()
max_or <- 50

df[["gba"]]$signif <- df[["gba"]]$FDR < 0.05
df[["gba"]]$OR_plot <- pmin(df[["gba"]]$OR, max_or)
df[["gba"]]$OR_inf  <- is.infinite(df[["gba"]]$OR)
df[["gba"]]$SYMBOL <- factor(df[["gba"]]$SYMBOL, levels = df[["gba"]]$SYMBOL[order(df[["gba"]]$OR)])
p[["gba"]] <- ggplot(df[["gba"]], aes(y = SYMBOL)) +
  geom_point(aes(x = OR_plot, color = signif)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = pmin(CI_upper, max_or), color = signif),
                 height = 0.2) +
  scale_color_manual(values = c("FALSE" = "black",
                                "TRUE"  = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  geom_text(
    data = subset(df[["gba"]], OR_inf),
    aes(x = max_or, label = "∞"),
    hjust = -0.2
  ) +
  coord_cartesian(xlim = c(0.2, max_or)) +
  theme_bw() +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL,
    title = NULL
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )

df[["sva"]]$signif <- df[["sva"]]$FDR < 0.05
df[["sva"]]$OR_plot <- pmin(df[["sva"]]$OR, max_or)
df[["sva"]]$OR_inf  <- is.infinite(df[["sva"]]$OR)
gene_order <- df[["gba"]]$SYMBOL[order(df[["gba"]]$OR)]
df[["sva"]]$SYMBOL <- factor(
  df[["sva"]]$SYMBOL,
  levels = gene_order
)
df[["sva"]]$label <- paste0(df[["sva"]]$SYMBOL, " - ", df[["sva"]]$Pos)

df[["sva"]]$label <- factor(
  df[["sva"]]$label,
  levels = df[["sva"]]$label[order(df[["sva"]]$SYMBOL)]
)

p[["sva"]] <- ggplot(df[["sva"]], aes(y = label)) +
  geom_point(aes(x = OR_plot, color = signif)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = pmin(CI_upper, max_or), color = signif),
                 height = 0.2) +
  scale_color_manual(values = c("FALSE" = "black",
                                "TRUE"  = "red")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  geom_text(
    data = subset(df[["sva"]], OR_inf),
    aes(x = max_or, label = "∞"),
    hjust = -0.2
  ) +
  coord_cartesian(xlim = c(0.2, max_or)) +
  theme_bw() +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL,
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )

pdf("figures/MCI_vs_nonMCI_OR_plot.pdf", width=25, height=15)
# png("figures/MCI_vs_nonMCI_OR_plot.png", width=25*300, height=15*300, res=300)
plot_grid(p[["gba"]], p[["sva"]],
        nrow=1, ncol=2, labels = c("A", "B"), rel_widths = c(0.8, 1.5))
dev.off()

#      SYMBOL No.of.variants.MCI
# 1     AGAP3                  2
# 2   ANKRD52                  1
# 3      ANO2                  1
# 4    ATP12A                  3
# 5  ATP6V1E2                  1
# 6    EIF4G1                  1
# 7     ENDOG                  2
# 8     EP400                  3
# 9     ERBB3                  1
# 10  GRAMD1A                  2
# 11      HK3                  5
# 12    IGFN1                  1
# 13 KIAA1328                  1
# 14  PIK3C2G                  4
# 15    PRKRA                  3
# 16      PRX                  1
# 17    REEP6                  1
# 18    SCN5A                  1
# 19     TBCC                  1
# 20    TMED3                  1
# 21 TNKS1BP1                  2
#                                                                                                                                                            Pos.MCI
# 1                                                                                                                            chr7:151141913:T>G;chr7:151143876:C>T
# 2                                                                                                                                               chr12:56248815:C>T
# 3                                                                                                                                                chr12:5744308:C>T
# 4                                                                                                         chr13:24690378:C>A;chr13:24698695:G>A;chr13:24706425:C>A
# 5                                                                                                                                                chr2:46512350:A>G
# 6                                                                                                                                               chr3:184315833:G>C
# 7                                                                                                                        chr9:128822601:CAGTA>C;chr9:128818882:G>C
# 8                                                                                                   chr12:132062583:AGC>A;chr12:132062586:AG>A;chr12:132053386:A>C
# 9                                                                                                                                               chr12:56101239:G>A
# 10                                                                                                                           chr19:35014388:G>A;chr19:35019400:C>T
# 11                                                                 chr5:176889700:G>GC;chr5:176881302:C>T;chr5:176887205:T>C;chr5:176889558:C>T;chr5:176891161:C>T
# 12                                                                                                                                              chr1:201213258:C>T
# 13                                                                                                                                              chr18:37067460:C>T
# 14                                                                             chr12:18338490:CCAAA>C;chr12:18562893:G>A;chr12:18594540:CAGTT>C;chr12:18399660:C>T
# 15 chr2:178441705:C>CAGTTTCCATAAATGACTCTAGCCTGCAAATTGTAGTATATTCTCTCTTA;chr2:178441705:C>CAGTTTCCATAAATGACTCTAGCCTGCAAATTGTAGTATATTCTCTCT;2:178444501-178447508:DEL
# 16                                                                                                                                              chr19:40397301:G>A
# 17                                                                                                                                               chr19:1496349:G>A
# 18                                                                                                                                               chr3:38633208:G>A
# 19                                                                                                                                               chr6:42745816:C>A
# 20                                                                                                                                              chr15:79313844:G>A
# 21                                                                                                                           chr11:57320343:C>T;chr11:57320509:C>G
#    AC.MCI AN.MCI RC.MCI gnomADg_AC_nfe.MCI gnomADg_AN_nfe.MCI
# 1       5    300    295                286             127954
# 2       3    150    147                 76              64574
# 3       2    150    148                 23              64574
# 4      15    450    435               1479             193598
# 5       4    150    146                158              64580
# 6       4    150    146                 11              62000
# 7       5    300    295               1006             128954
# 8       4    450    446                 15             192132
# 9       2    150    148                  2              64544
# 10      6    300    294               1914             129136
# 11     12    750    738               3833             322812
# 12      2    150    148                 16              64548
# 13      2    150    148                 22              64572
# 14      6    600    594                560             258146
# 15     18    450    432                 45             178962
# 16      2    150    148                 14              64554
# 17      2    150    148                 32              64566
# 18      3    150    147                 81              64562
# 19      3    150    147                 99              64582
# 20      2    150    148                 36              64584
# 21     12    300    288               1479             129146
#    gnomADg_RC_nfe.MCI gnomADg_AC.MCI gnomADg_AN.MCI gnomADg_RC.MCI
# 1              127668            541         283946         283405
# 2               64498            114         143290         143176
# 3               64551            618         143302         142684
# 4              192119           3596         429122         425526
# 5               64422            284         143214         142930
# 6               61989             37         138228         138191
# 7              127948           1559         285962         284403
# 8              192117             20         424654         424634
# 9               64542             60         143162         143102
# 10             127222           2851         286562         283711
# 11             318979           6356         716434         710078
# 12              64532             56         143196         143140
# 13              64550           4644         143176         138532
# 14             257586            772         572700         571928
# 15             178917             82         394496         394414
# 16              64540           2898         143232         140334
# 17              64534             99         143318         143219
# 18              64481           4326         143218         138892
# 19              64483            170         143364         143194
# 20              64548            977         143306         142329
# 21             127667           2709         286610         283901
#                                                                                                                                                                                                           Pos.nonMCI
# 1                                                                                                                                                                              chr7:151141913:T>G;chr7:151143876:C>T
# 2                                                                                                                                                                                                 chr12:56248815:C>T
# 3                                                                                                                                                                                                  chr12:5744308:C>T
# 4                                                                                                                                                           chr13:24690378:C>A;chr13:24698695:G>A;chr13:24706425:C>A
# 5                                                                                                                                                                                                  chr2:46512350:A>G
# 6                                                                                                                                                                                                 chr3:184315833:G>C
# 7                                                                                                                                                                          chr9:128822601:CAGTA>C;chr9:128818882:G>C
# 8                                                                                                                                                                                                chr12:132053386:A>C
# 9                                                                                                                                                                                                 chr12:56101239:G>A
# 10                                                                                                                                                                             chr19:35014388:G>A;chr19:35019400:C>T
# 11                                                                                                                   chr5:176889700:G>GC;chr5:176881302:C>T;chr5:176887205:T>C;chr5:176889558:C>T;chr5:176891161:C>T
# 12                                                                                                                                                                            chr1:201213258:C>T;chr1:201215703:GC>G
# 13                                                                                                                                                                                                chr18:37067460:C>T
# 14                                                                                                                                                      chr12:18562893:G>A;chr12:18594540:CAGTT>C;chr12:18399660:C>T
# 15 chr2:178432255:C>CTATATCCAAATA;chr2:178441705:C>CAGTTTCCATAAATGACTCTAGCCTGCAAATTGTAGTATATTCTCTCTTA;chr2:178441705:C>CAGTTTCCATAAATGACTCTAGCCTGCAAATTGTAGTATATTCTCTCT;chr2:178450999:G>A;2:178444501-178447508:DEL
# 16                                                                                                                                                                                                chr19:40397301:G>A
# 17                                                                                                                                                                                                 chr19:1496349:G>A
# 18                                                                                                                                                                                                 chr3:38633208:G>A
# 19                                                                                                                                                                                                 chr6:42745816:C>A
# 20                                                                                                                                                                                                chr15:79313844:G>A
# 21                                                                                                                                                                             chr11:57320343:C>T;chr11:57320509:C>G
#    AC.nonMCI AN.nonMCI RC.nonMCI gnomADg_AC_nfe.nonMCI gnomADg_AN_nfe.nonMCI
# 1          8       804       796                   286                127954
# 2          0       150       150                    76                 64574
# 3          0       150       150                    23                 64574
# 4         14      1206      1192                  1479                193598
# 5          3       402       399                   158                 64580
# 6         18       402       384                    11                 62000
# 7          9       804       795                  1006                128954
# 8          7       402       395                     1                 63116
# 9          0       150       150                     2                 64544
# 10         8       804       796                  1914                129136
# 11        26      2010      1984                  3833                322812
# 12         4       804       800                    49                129108
# 13         0       150       150                    22                 64572
# 14         6      1206      1200                   549                193638
# 15        19      2010      1991                   383                302660
# 16         1       402       401                    14                 64554
# 17         1       402       401                    32                 64566
# 18         1       402       401                    81                 64562
# 19         3       402       399                    99                 64582
# 20         0       150       150                    36                 64584
# 21        14       804       790                  1479                129146
#    gnomADg_RC_nfe.nonMCI gnomADg_AC.nonMCI gnomADg_AN.nonMCI gnomADg_RC.nonMCI
# 1                 127668               541            283946            283405
# 2                  64498               114            143290            143176
# 3                  64551               618            143302            142684
# 4                 192119              3596            429122            425526
# 5                  64422               284            143214            142930
# 6                  61989                37            138228            138191
# 7                 127948              1559            285962            284403
# 8                  63115                 4            138396            138392
# 9                  64542                60            143162            143102
# 10                127222              2851            286562            283711
# 11                318979              6356            716434            710078
# 12                129059                96            286466            286370
# 13                 64550              4644            143176            138532
# 14                193089               758            429674            428916
# 15                302277               866            669848            668982
# 16                 64540              2898            143232            140334
# 17                 64534                99            143318            143219
# 18                 64481              4326            143218            138892
# 19                 64483               170            143364            143194
# 20                 64548               977            143306            142329
# 21                127667              2709            286610            283901
#    p_value_NFE_gnomADg.MCI FDR_NFE_gnomADg.MCI OR_NFE_gnomADg.MCI
# 1             8.626756e-06        1.811619e-04           7.565959
# 2             1.603758e-06        3.367891e-05          17.319549
# 3             9.525519e-07        2.000359e-05          37.926557
# 4             1.333238e-08        2.799799e-07           4.479238
# 5             2.550236e-06        5.355496e-05          11.170799
# 6             1.279486e-17        2.686920e-16         154.393524
# 7             8.932267e-02        1.000000e+00           2.155676
# 8             4.458241e-17        9.362306e-16         114.868161
# 9             1.383541e-09        2.905437e-08         436.094595
# 10            4.603706e-01        1.000000e+00           1.356514
# 11            2.994319e-01        1.000000e+00           1.353157
# 12            1.163034e-07        2.442371e-06          54.503378
# 13            7.337365e-07        1.540847e-05          39.649877
# 14            1.961348e-04        4.118831e-03           4.646212
# 15            7.158814e-73        1.503351e-71         165.663889
# 16            5.513276e-08        1.157788e-06          62.297297
# 17            6.606300e-06        1.387323e-04          27.252534
# 18            2.658049e-06        5.581903e-05          16.246158
# 19            1.233295e-05        2.589919e-04          13.292723
# 20            1.301790e-05        2.733759e-04          24.229730
# 21            1.508341e-05        3.167515e-04           3.596659
#    CI_lower_NFE_gnomADg.MCI CI_upper_NFE_gnomADg.MCI EF_NFE_gnomADg.MCI
# 1                 3.1022914                18.452084          0.8678291
# 2                 5.4022471                55.526296          0.9422618
# 3                 8.8617626               162.318016          0.9736333
# 4                 2.6703123                 7.513568          0.7767477
# 5                 4.0868541                30.533695          0.9104809
# 6                48.6010638               490.469929          0.9935230
# 7                 0.8886910                 5.228971          0.5361084
# 8                37.9748680               347.458601          0.9912944
# 9                61.0219462              3116.558998          0.9977069
# 10                0.6037279                 3.047945          0.2628162
# 11                0.7642710                 2.395791          0.2609874
# 12               12.4216397               239.148642          0.9816525
# 13                9.2402908               170.136718          0.9747792
# 14                2.0700887                10.428194          0.7847709
# 15               95.1304565               288.493560          0.9939637
# 16               14.0348497               276.522609          0.9839479
# 17                6.4720085               114.755813          0.9633062
# 18                5.0742334                52.015274          0.9384470
# 19                4.1672102                42.401626          0.9247709
# 20                5.7809440               101.554314          0.9587284
# 21                2.0142848                 6.422108          0.7219642
#    p_value_NFE_gnomADg.nonMCI FDR_NFE_gnomADg.nonMCI OR_NFE_gnomADg.nonMCI
# 1                3.086733e-05           4.938772e-04             4.4863478
# 2                         NaN                    NaN             0.0000000
# 3                         NaN                    NaN             0.0000000
# 4                1.178203e-01           1.000000e+00             1.5256465
# 5                5.547985e-02           8.876777e-01             3.0656705
# 6                2.787044e-47           4.459270e-46           264.1576705
# 7                2.789773e-01           1.000000e+00             1.4398290
# 8                5.414889e-11           8.663823e-10          1118.4936709
# 9                         NaN                    NaN             0.0000000
# 10               2.572251e-01           1.000000e+00             0.6680319
# 11               6.615641e-01           1.000000e+00             1.0905735
# 12               7.578835e-07           1.212614e-05            13.1692857
# 13                        NaN                    NaN             0.0000000
# 14               1.701225e-01           1.000000e+00             1.7585519
# 15               1.214419e-17           1.943070e-16             7.5316247
# 16               1.845009e-02           2.952014e-01            11.4962594
# 17               1.121379e-01           1.000000e+00             5.0291459
# 18               4.960751e-01           1.000000e+00             1.9851913
# 19               6.912681e-03           1.106029e-01             4.8973191
# 20                        NaN                    NaN             0.0000000
# 21               1.165913e-01           1.000000e+00             1.5297182
#    CI_lower_NFE_gnomADg.nonMCI CI_upper_NFE_gnomADg.nonMCI
# 1                    2.2144451                    9.089102
# 2                    0.0000000                         NaN
# 3                    0.0000000                         NaN
# 4                    0.8985626                    2.590356
# 5                    0.9740765                    9.648457
# 6                  123.9353893                  563.029456
# 7                    0.7442196                    2.785612
# 8                  137.2873242                 9112.480699
# 9                    0.0000000                         NaN
# 10                   0.3324324                    1.342428
# 11                   0.7397056                    1.607870
# 12                   4.7412666                   36.578852
# 13                   0.0000000                         NaN
# 14                   0.7850279                    3.939357
# 15                   4.7414094                   11.963820
# 16                   1.5081062                   87.635727
# 17                   0.6855325                   36.894395
# 18                   0.2756044                   14.299427
# 19                   1.5462837                   15.510565
# 20                   0.0000000                         NaN
# 21                   0.8995596                    2.601315
#    EF_NFE_gnomADg.nonMCI      p_value          FDR        OR  CI_lower
# 1             0.77710154 3.626890e-01 1.0000000000 1.6864407 0.5473263
# 2                   -Inf          NaN          NaN       Inf       NaN
# 3                   -Inf          NaN          NaN       Inf       NaN
# 4             0.34454018 4.158209e-03 0.0665313385 2.9359606 1.4055764
# 5             0.67380708 9.304146e-02 1.0000000000 3.6438356 0.8058203
# 6             0.99621438 3.386402e-01 1.0000000000 0.5844749 0.1945425
# 7             0.30547306 4.726272e-01 1.0000000000 1.4971751 0.4976833
# 8             0.99910594 2.801261e-01 1.0000000000 0.5060858 0.1470496
# 9                   -Inf          NaN          NaN       Inf       NaN
# 10           -0.49693449 1.931772e-01 1.0000000000 2.0306122 0.6986524
# 11            0.08305128 5.395409e-01 1.0000000000 1.2407755 0.6228268
# 12            0.92406574 2.534606e-01 1.0000000000 2.7027027 0.4905667
# 13                  -Inf          NaN          NaN       Inf       NaN
# 14            0.43135031 2.249763e-01 1.0000000000 2.0202020 0.6487726
# 15            0.86722652 9.696763e-06 0.0001551482 4.3662281 2.2725000
# 16            0.91301518 1.689580e-01 1.0000000000 5.4189189 0.4877296
# 17            0.80115908 1.689580e-01 1.0000000000 5.4189189 0.4877296
# 18            0.49627022 6.964786e-02 1.0000000000 8.1836735 0.8445527
# 19            0.79580665 2.245595e-01 1.0000000000 2.7142857 0.5417564
# 20                  -Inf          NaN          NaN       Inf       NaN
# 21            0.34628481 3.230210e-02 0.5168335518 2.3511905 1.0748143
#     CI_upper         EF       p_diff     FDR_diff IS.SIG
# 1   5.196320  0.4070352 3.677328e-01 4.525942e-01     ns
# 2        Inf        NaN          NaN          NaN     ns
# 3        Inf        NaN          NaN          NaN     ns
# 4   6.132619  0.6593960 4.342090e-03 3.473672e-02     ns
# 5  16.477045  0.7255639 9.653848e-02 2.206594e-01     ns
# 6   1.755970 -0.7109375 4.461270e-01 5.056284e-01     ns
# 7   4.503935  0.3320755 4.740266e-01 5.056284e-01     ns
# 8   1.741745 -0.9759494 5.999743e-02 2.206594e-01     ns
# 9        Inf        NaN          NaN          NaN     ns
# 10  5.901914  0.5075377 1.939732e-01 2.821428e-01     ns
# 11  2.471833  0.1940524 5.404043e-01 5.404043e-01     ns
# 12 14.890129  0.6300000 1.214077e-01 2.428154e-01     ns
# 13       Inf        NaN          NaN          NaN     ns
# 14  6.290673  0.5050000 9.541039e-02 2.206594e-01     ns
# 15  8.388976  0.7709694 5.030296e-17 8.048473e-16     ns
# 16 60.206894  0.8154613 1.886004e-01 2.821428e-01     ns
# 17 60.206894  0.8154613 1.776843e-01 2.821428e-01     ns
# 18 79.299391  0.8778055 7.222164e-02 2.206594e-01     ns
# 19 13.599003  0.6315789 2.314175e-01 3.085567e-01     ns
# 20       Inf        NaN          NaN          NaN     ns
# 21  5.143304  0.5746835 3.304540e-02 1.762421e-01     ns

### This script to calculate OR and adjusted P-value for structural variants in MCI vs. external databases
df <- data.frame(
    SYMBOL = c("PRKRA", "PIK3C2G", "ATP6V1E2", "TNKS1BP1"),
    Pos = c("chr2:178444501-178447508:DEL", "chr12:18570336-18572160:DEL", "chr2:46524048-46524172:DEL", "chr11:57304632-57304828:DEL"),
    AC = c(14, 1, 1, 1),
    AN = rep((n.MCI * 2), 4),
    gnomADg_AC_nfe = c(10, 753, 5, 650),
    gnomADg_AN_nfe = c(49918, 59084, 58892, 58178),
    gnomADg_AC = c(37, 4465, 484, 888),
    gnomADg_AN = c(108658, 126064, 125610, 123942)
)

df$RC <- df$AN - df$AC
df$gnomADg_RC_nfe <- df$gnomADg_AN_nfe
df$gnomADg_RC <- df$gnomADg_AN - df$gnomADg_AC

### NFE gnomeADg
t.gnomADg.nfe <- apply(df[, c("AC", "RC", "gnomADg_AC_nfe", "gnomADg_RC_nfe")], 1, function(x) {
    t_res <- matrix(x, nrow = 2, byrow=T)
    OR <- (x[1] * x[4]) / (x[2] * x[3])
    logOR <- log(OR)
    SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
    z <- logOR / SE_logOR
    p_value <- 2 * pnorm(-abs(z))
    CI_lower <- exp(logOR - 1.96 * SE_logOR)
    CI_upper <- exp(logOR + 1.96 * SE_logOR)
    return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
})
t.gnomADg.nfe <- as.data.frame(t(t.gnomADg.nfe))
adj.pVal.gnomADg.nfe <- p.adjust(t.gnomADg.nfe$pVal, method = "bonferroni")
EF.gnomADg.nfe <- (t.gnomADg.nfe$OR-1)/t.gnomADg.nfe$OR

### ALL gnomeADg
t.gnomADg <- apply(df[, c("AC", "RC", "gnomADg_AC", "gnomADg_RC")], 1, function(x) {
    t_res <- matrix(x, nrow = 2, byrow=T)
    OR <- (x[1] * x[4]) / (x[2] * x[3])
    logOR <- log(OR)
    SE_logOR <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
    z <- logOR / SE_logOR
    p_value <- 2 * pnorm(-abs(z))
    CI_lower <- exp(logOR - 1.96 * SE_logOR)
    CI_upper <- exp(logOR + 1.96 * SE_logOR)
    return(c(pVal=p_value, OR=OR, CI_lower=CI_lower, CI_upper=CI_upper))
})
t.gnomADg <- as.data.frame(t(t.gnomADg))
adj.pVal.gnomADg <- p.adjust(t.gnomADg$pVal, method = "bonferroni")
EF.gnomADg <- (t.gnomADg$OR-1)/t.gnomADg$OR

df.res <- cbind(df, 
                p_value_NFE_gnomADg=t.gnomADg.nfe$pVal, 
                FDR_NFE_gnomADg=adj.pVal.gnomADg.nfe, 
                OR_NFE_gnomADg=t.gnomADg.nfe$OR,
                CI_lower_NFE_gnomADg=t.gnomADg.nfe$CI_lower,
                CI_upper_NFE_gnomADg=t.gnomADg.nfe$CI_upper,
                EF_NFE_gnomADg=EF.gnomADg.nfe,   
                p_value_ALL_gnomADg=t.gnomADg$pVal, 
                FDR_ALL_gnomADg=adj.pVal.gnomADg,
                OR_ALL_gnomADg=t.gnomADg$OR,
                CI_lower_ALL_gnomADg=t.gnomADg$CI_lower,
                CI_upper_ALL_gnomADg=t.gnomADg$CI_upper,
                EF_ALL_gnomADg=EF.gnomADg,
                IS.SIG=rep("ns", nrow(df)))

df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg > 0.05] <- "NFE"
df.res$IS.SIG[df.res$OR_ALL_gnomADg > 1 & df.res$FDR_ALL_gnomADg < 0.05 & df.res$FDR_NFE_gnomADg > 0.05] <- "ALL"
df.res$IS.SIG[df.res$OR_NFE_gnomADg > 1 & df.res$OR_ALL_gnomADg > 1 & df.res$FDR_NFE_gnomADg < 0.05 & df.res$FDR_ALL_gnomADg < 0.05] <- "NFE/ALL"

df.res <- df.res[order(df.res$FDR_NFE_gnomADg), ]
write.csv(df.res, "GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.SV.GBA.csv")