
library(synaptome.db) # nolint: object_name_linter.
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/")

LoF_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.LoF.GBA.v3.csv", header=T)
missense_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.missense.GBA.v3.csv", header=T)
LoF_genes <- LoF_dt$SYMBOL
missense_genes <- missense_dt$SYMBOL

# Extract synaptic genes
LoF_synap_dt <- as.data.frame(getGeneDiseaseByName(LoF_genes))
missense_synap_dt <- as.data.frame(getGeneDiseaseByName(missense_genes))

brain_disease <- unique(c(LoF_synap_dt$Description, missense_synap_dt$Description))
# target_brain_disease <- c("brain_disease", "central_nervous_system_disease", "nervous_system_disease", "schizophrenia", "psychotic_disorder","cognitive_disorder",
#     "disease_of_mental_health", "neurodegenerative_disease", "Parkinson's_disease", "dystonia", "neuropathy", "nervous_system_cancer", "Lewy_body_dementia",
#     "dementia", "neuromuscular_disease", "autonomic_nervous_system_neoplasm", "neuroblastoma", "peripheral_nervous_system_neoplasm", "sensory_peripheral_neuropathy")
target_brain_disease <- c("brain_disease", "nervous_system_disease", "central_nervous_system_disease", 
    "neurodegenerative_disease", "dementia", "Parkinson's_disease", "Lewy_body_dementia",
    "disease_of_mental_health", "cognitive_disorder", "schizophrenia", "psychotic_disorder",
    "neuropathy", "sensory_peripheral_neuropathy",
    "nervous_system_cancer", "autonomic_nervous_system_neoplasm", "neuroblastoma", "peripheral_nervous_system_neoplasm",
    "neuromuscular_disease", "dystonia")

group_colors <- c(
  "Broad umbrella terms" = "#1f78b4",
  "Neurodegenerative diseases" = "#e31a1c",
  "Cognitive & Psychiatric disorders" = "#ff7f00",
  "Peripheral & sensory nerve disorders" = "#6a3d9a",
  "Nervous system cancers" = "#b15928",
  "Neuromuscular disorders" = "#a6cee3",
  "Other" = "grey60"
)

BD_dt <- rbind(LoF_synap_dt[LoF_synap_dt$Description %in% target_brain_disease, ], missense_synap_dt[missense_synap_dt$Description %in% target_brain_disease, ])                               

dot_data <- BD_dt %>%
  select(HumanName, Description) %>%
  rename(gene = HumanName, disease = Description)

dot_data <- dot_data %>%
  mutate(disease_group = case_when(
    disease %in% c("brain_disease", "nervous_system_disease", "central_nervous_system_disease") ~ "Broad umbrella terms",
    disease %in% c("neurodegenerative_disease", "dementia", "Parkinson's_disease", "Lewy_body_dementia") ~ "Neurodegenerative diseases",
    disease %in% c("disease_of_mental_health", "cognitive_disorder", "schizophrenia", "psychotic_disorder") ~ "Cognitive & Psychiatric disorders",
    disease %in% c("neuropathy", "sensory_peripheral_neuropathy") ~ "Peripheral & sensory nerve disorders",
    disease %in% c("nervous_system_cancer", "autonomic_nervous_system_neoplasm", "neuroblastoma", "peripheral_nervous_system_neoplasm") ~ "Nervous system cancers",
    disease %in% c("neuromuscular_disease", "dystonia") ~ "Neuromuscular disorders",
    TRUE ~ "Other"
  ))

# Sort genes by number of associated diseases
gene_order <- dot_data %>%
    count(gene) %>%
    arrange(desc(n)) %>%
    pull(gene)

dot_data$gene <- factor(dot_data$gene, levels = gene_order)
dot_data$disease <- factor(dot_data$disease, levels = target_brain_disease)
dot_data$disease_group <- factor(dot_data$disease_group, levels = c(
    "Broad umbrella terms", "Neurodegenerative diseases", "Cognitive & Psychiatric disorders",
    "Peripheral & sensory nerve disorders", "Nervous system cancers", "Neuromuscular disorders", "Other"
))

# Extract sample data for bar plot
LoF_samples <- LoF_dt[, c("SYMBOL", "samples_het", "samples_hom")]
missense_samples <- missense_dt[, c("SYMBOL", "samples_het", "samples_hom")]
sample_data <- rbind(LoF_samples, missense_samples)
sample_data <- sample_data[sample_data$SYMBOL %in% dot_data$gene, ]
sample_data$SYMBOL <- factor(sample_data$SYMBOL, levels = gene_order)
sample_data <- sample_data[order(sample_data$SYMBOL), ]

#     SYMBOL              samples_het samples_hom
# 21   ERBB3                AT07|AT13        <NA>
# 16  EIF4G1      AT13|AT21|LV45|UG23        <NA>
# 1    PRKRA      BL13|LV29;RG11|RG42       NA;NA
# 40   SCN5A           LV51|UG36|UG72        <NA>
# 11   ENDOG           LV34|UG17|UG23        <NA>
# 5  PIK3C2G AT14;BL13|RG12|RG27;RG56    NA;NA;NA
# 27     PRX                AT16|UG51        <NA>
# 6  GRAMD1A      BL13|BL20|BL29|UG28        <NA>
# 48    TBCC           RG27|UG10|UG65        <NA>
# 33    ANO2                BL33|UG51        <NA>
# 3    IGFN1                UG69|UG72        <NA>
# 37 ANKRD52           BL20|UG01|UG65        <NA>

total_samples <- c()
for (row in 1:nrow(sample_data)) {
    s1 <- unique(unlist(strsplit(sample_data$samples_het[row], "[|;]")))
    s2 <- unique(unlist(strsplit(sample_data$samples_hom[row], "[|;]")))
    s1 <- s1[!is.na(s1) & s1 != "NA"]
    s2 <- s2[!is.na(s2) & s2 != "NA"]
    total_samples <- unique(c(total_samples, s1, s2))
    sample_data$cum_freq[row] <- length(total_samples)/75
}

dot_plot <- ggplot(dot_data, aes(x = gene, y = disease)) +
    geom_point(aes(color = disease_group), size = 6) +
    scale_color_manual(values = group_colors, name = "Associated Diseases (HDO)") +
    theme_bw() +
    theme(axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = group_colors[dot_data$disease_group[match(levels(dot_data$disease),
                                                        dot_data$disease)]], size = 12),
        legend.position = "right",
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))

bar_data <- data.frame(
  gene = unique(dot_data$gene),
  n_samples = c(4, 4, 3, 5, 2, 4, 3, 2, 3, 2, 3, 2)
)
#  [1] PRKRA   GRAMD1A ENDOG   PIK3C2G IGFN1   EIF4G1  SCN5A   ERBB3   TBCC   
# [10] PRX     ANKRD52 ANO2   

bar_data <- merge(sample_data, bar_data, by.x="SYMBOL", by.y="gene")
bar_data$SYMBOL <- factor(bar_data$SYMBOL, levels = gene_order)
bar_data <- bar_data[order(bar_data$SYMBOL),]
max_cum <- max(bar_data$cum_freq, na.rm = TRUE)
scale_factor <- max(bar_data$n_samples) / max_cum

bar_plot <- ggplot(bar_data, aes(x = SYMBOL)) +
  geom_col(aes(y = n_samples), fill = "grey", width=0.6) +
  geom_line(aes(y = cum_freq * scale_factor, group = 1), color = "skyblue", size = 1) +
  geom_point(aes(y = cum_freq * scale_factor ), color = "skyblue", size=6) +
  scale_y_continuous(
    name = "# samples \nwith variants in gene",
    limits = c(0,6),
    sec.axis = sec_axis(~ . *100 / scale_factor, name = "Cumulative (%)")
  ) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 14),
            axis.title.y = element_text(size=16))

dir.create("figures", showWarnings = FALSE)
pdf("figures/synaptome.pdf", height = 8, width = 15)
png("figures/synaptome.png", height = 8*300, width = 15*300, res = 300)
plot_grid(bar_plot, dot_plot, ncol = 1, align = "v", axis = "lr", rel_heights = c(1,2))
dev.off()

# Synaptic genes associated with THI, GUF, MOCA, PHQ9 levels
synap_gba <- read.table("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.GBA_genes.v3.txt", header=TRUE)
BD_gba_genes <- unique(BD_dt$HumanName)          

gba_genes <- list()
# THI
for (q in c("mild", paste0("g", c(3:5)))) {
  gene <- c()
  for (c in c("LoF", "missense")) {
    dt <- read.csv(paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.csv"), header=T)
    gene <- c(gene, dt$SYMBOL)
  }
  gba_genes[["THI"]][[q]] <- unique(gene)
}
# GUF
for (q in paste0("g", c(1:4))) {
  gene <- c()
  for (c in c("LoF", "missense")) {
    dt <- read.csv(paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.csv"), header=T)
    gene <- c(gene, dt$SYMBOL)
  }
  gba_genes[["GUEF"]][[q]] <- unique(gene)
}
# PHQ9
for (q in paste0("g", c(2:4))) {
  gene <- c()
  for (c in c("LoF", "missense")) {
    dt <- read.csv(paste0("GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_", q, ".", c, ".GBA.csv"), header=T)
    gene <- c(gene, dt$SYMBOL)
  }
  gba_genes[["PHQ9"]][[q]] <- unique(gene)
}
# # MOCA
# for (q in "normal")) {
#   gene <- c()
#   for (c in c("LoF", "missense")) {
#     dt <- read.csv(paste0("GBA_MOCA/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_", q, ".", c, ".GBA.csv"), header=T)
#     gene <- c(gene, dt$SYMBOL)
#   }
#   gba_genes[["MOCA"]][[q]] <- unique(gene)
# }

synap_gba_dt <- data.frame(matrix(nrow=0, ncol=2))
colnames(synap_gba_dt) <- c("Gene", "Feature")
for (g in synap_gba$GBA_gene) {
  for (t in c("THI", "GUEF", "PHQ9", "MOCA")) {
    for (q in names(gba_genes[[t]])) {
      if (g %in% gba_genes[[t]][[q]]) {
        synap_gba_dt <- rbind(synap_gba_dt, data.frame(Gene=g, Feature=paste0(t, "_", q)))
      }
    }
  }
}


gba_dt <- rbind(LoF_dt, missense_dt)
gba_dt <- gba_dt[gba_dt$SYMBOL %in% synap_gba_dt$Gene, ]
gba_dt$gene_type <- ifelse(gba_dt$SYMBOL %in% BD_gba_genes, "BD-associated", "Non-BD-associated")
gba_dt$gene_type <- factor(gba_dt$gene_type, levels = c("BD-associated", "Non-BD-associated"))
gene_order <- gba_dt[order(gba_dt$gene_type, gba_dt$FDR_NFE_gnomADg), "SYMBOL"]
gene_order <- unique(gene_order)

synap_gba_dt$Gene <- factor(synap_gba_dt$Gene, levels = gene_order)
synap_gba_dt$Feature <- factor(synap_gba_dt$Feature, levels = c(
  "THI_mild", "THI_g3", "THI_g4", "THI_g5",
  "GUEF_g1", "GUEF_g2", "GUEF_g3", "GUEF_g4",
  "PHQ9_g2", "PHQ9_g3", "PHQ9_g4"
))
synap_gba_dt$Gene_type <- ifelse(synap_gba_dt$Gene %in% BD_gba_genes, "BD-associated", "Non-BD-associated")
synap_gba_dt$Gene_type <- factor(synap_gba_dt$Gene_type, levels = c("BD-associated", "Non-BD-associated"))
synap_gba_dt <- synap_gba_dt %>%
  mutate(Feature_group = case_when(
    Feature %in% c("THI_mild", "THI_g3", "THI_g4", "THI_g5") ~ "THI",
    Feature %in% c("GUEF_g1", "GUEF_g2", "GUEF_g3", "GUEF_g4") ~ "GUEF",
    Feature %in% c("PHQ9_g2", "PHQ9_g3", "PHQ9_g4") ~ "PHQ9",
    TRUE ~ "Other"
  ))

group_colors <- c(
  THI = "#4D4D4D",
  GUEF = "#969696",
  PHQ9 = "#D9D9D9"
)

pdf("figures/feature_THI_GUF_PHQ9.pdf", height = 6, width = 15)
ggplot(synap_gba_dt, aes(x = Gene, y = Feature)) +
    geom_point(aes(color = Feature_group), size = 6) +
    scale_y_discrete(drop = FALSE) +        
    scale_color_manual(values = group_colors) +
    theme_bw() +
    theme(axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.position = "right",
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))
dev.off()

# Non-synaptic genes associated with THI, GUF, MOCA, PHQ9 levels
gba_dt <- rbind(LoF_dt, missense_dt)
all_gba_genes <- unique(gba_dt$SYMBOL)
non_synap_gba <- all_gba_genes[!(all_gba_genes %in% synap_gba$GBA_gene)]

non_synap_gba_dt <- data.frame(matrix(nrow=0, ncol=2))
colnames(non_synap_gba_dt) <- c("Gene", "Feature")
for (g in non_synap_gba) {
  for (t in c("THI", "GUEF", "PHQ9", "MOCA")) {
    for (q in names(gba_genes[[t]])) {
      if (g %in% gba_genes[[t]][[q]]) {
        non_synap_gba_dt <- rbind(non_synap_gba_dt, data.frame(Gene=g, Feature=paste0(t, "_", q)))
      }
    }
  }
}


gba_dt <- rbind(LoF_dt, missense_dt)
gba_dt <- gba_dt[gba_dt$SYMBOL %in% non_synap_gba_dt$Gene, ]
gene_order <- gba_dt[order(gba_dt$FDR_NFE_gnomADg), "SYMBOL"]
gene_order <- unique(gene_order)

non_synap_gba_dt$Gene <- factor(non_synap_gba_dt$Gene, levels = gene_order)
non_synap_gba_dt$Feature <- factor(non_synap_gba_dt$Feature, levels = c(
  "THI_mild", "THI_g3", "THI_g4", "THI_g5",
  "GUEF_g1", "GUEF_g2", "GUEF_g3", "GUEF_g4",
  "PHQ9_g2", "PHQ9_g3", "PHQ9_g4"
))
non_synap_gba_dt <- non_synap_gba_dt %>%
  mutate(Feature_group = case_when(
    Feature %in% c("THI_mild", "THI_g3", "THI_g4", "THI_g5") ~ "THI",
    Feature %in% c("GUEF_g1", "GUEF_g2", "GUEF_g3", "GUEF_g4") ~ "GUEF",
    Feature %in% c("PHQ9_g2", "PHQ9_g3", "PHQ9_g4") ~ "PHQ9",
    TRUE ~ "Other"
  ))

group_colors <- c(
  THI = "#4D4D4D",
  GUEF = "#969696",
  PHQ9 = "#D9D9D9"
)

pdf("figures/feature_THI_GUF_PHQ9_non_synaptic.pdf", height = 6, width = 15)
ggplot(non_synap_gba_dt, aes(x = Gene, y = Feature)) +
    geom_point(aes(color = Feature_group), size = 6) +
    scale_y_discrete(drop = FALSE) +        
    scale_color_manual(values = group_colors) +
    theme_bw() +
    theme(axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.position = "right",
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))
dev.off()

# THI vs GUF
gba_genes <- list()
# THI
for (q in c("mild", paste0("g", c(3:5)))) {
  gene <- c()
  for (c in c("LoF", "missense")) {
    dt <- read.csv(paste0("GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_", q, ".", c, ".GBA.csv"), header=T)
    gene <- c(gene, dt$SYMBOL)
  }
  gba_genes[["THI"]][[q]] <- unique(gene)
}
# GUF
for (q in paste0("g", c(1:4))) {
  gene <- c()
  for (c in c("LoF", "missense")) {
    dt <- read.csv(paste0("GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_", q, ".", c, ".GBA.csv"), header=T)
    gene <- c(gene, dt$SYMBOL)
  }
  gba_genes[["GUEF"]][[q]] <- unique(gene)
}



pl1_i <- intersect(unique(c(gba_genes[["THI"]][["mild"]], gba_genes[["THI"]][["g3"]])), 
          unique(c(gba_genes[["GUEF"]][["g1"]], gba_genes[["GUEF"]][["g2"]])))
pl2_i <- intersect(unique(c(gba_genes[["THI"]][["g4"]], gba_genes[["THI"]][["g5"]])), 
          unique(c(gba_genes[["GUEF"]][["g3"]], gba_genes[["GUEF"]][["g4"]])))
pl3_i <- intersect(unique(c(gba_genes[["THI"]][["g4"]], gba_genes[["THI"]][["g5"]])), 
          unique(c(gba_genes[["GUEF"]][["g1"]], gba_genes[["GUEF"]][["g2"]])))
pl4_i <- intersect(unique(c(gba_genes[["THI"]][["mild"]], gba_genes[["THI"]][["g3"]])), 
          unique(c(gba_genes[["GUEF"]][["g3"]], gba_genes[["GUEF"]][["g4"]])))

library(VennDiagram)
pl1 <- list("Mild-moderate THI" = unique(c(gba_genes[["THI"]][["mild"]], gba_genes[["THI"]][["g3"]])),
                "Mild-moderate GUF" = unique(c(gba_genes[["GUEF"]][["g1"]], gba_genes[["GUEF"]][["g2"]])))

pl2 <- list("Severe THI" = unique(c(gba_genes[["THI"]][["g4"]], gba_genes[["THI"]][["g5"]])),
                "Severe GUF" = unique(c(gba_genes[["GUEF"]][["g3"]], gba_genes[["GUEF"]][["g4"]])))

pl3 <- list("Severe THI" = unique(c(gba_genes[["THI"]][["g4"]], gba_genes[["THI"]][["g5"]])),
                "Mild-moderate GUF" = unique(c(gba_genes[["GUEF"]][["g1"]], gba_genes[["GUEF"]][["g2"]])))

pl4 <- list("Mild-moderate THI" = unique(c(gba_genes[["THI"]][["mild"]], gba_genes[["THI"]][["g3"]])),
                "Severe GUF" = unique(c(gba_genes[["GUEF"]][["g3"]], gba_genes[["GUEF"]][["g4"]])))

p1 <- venn.diagram(pl1, filename = NULL, disable.logging = TRUE, 
             main=NULL, print.mode="raw",
             fill = c("darkcyan", "grey"), lwd=1, cex=1.5, cat.cex=1.2)
p2 <- venn.diagram(pl2, filename = NULL, disable.logging = TRUE, 
             main=NULL, print.mode="raw",
             fill = c("darkcyan", "grey"), lwd=1, cex=1.5, cat.cex=1.2)
p3 <- venn.diagram(pl3, filename = NULL, disable.logging = TRUE, 
             main=NULL, print.mode="raw",
             fill = c("darkcyan", "grey"), lwd=1, cex=1.5, cat.cex=1.2)
p4 <- venn.diagram(pl4, filename = NULL, disable.logging = TRUE, 
             main=NULL, print.mode="raw",
             fill = c("darkcyan", "grey"), lwd=1, cex=1.5, cat.cex=1.2)

# png("figures/THI_vs_GUF.png", height = 10*300, width = 10*300, res = 300)
pdf("figures/THI_vs_GUF.pdf", height = 10, width = 10)
plot_grid(p1, p2, p3, p4, ncol=2, nrow=2)
dev.off()