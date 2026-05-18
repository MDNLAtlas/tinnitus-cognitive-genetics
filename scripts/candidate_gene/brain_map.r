setwd("/Users/phammaitam/Documents/research_projects/MD/Tinnitus-MCD/output/")

library(ggplot2)
library(ggseg)
library(gridExtra)
library(dplyr)

# Gene expression
dt_1 <- tibble(
  region = rep(c("entorhinal", 
             "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
             "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medical orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis", "precentral gyrus",
             "parahippocampal"), 2),
  log2FC = c(-0.061,
             rep(0.483, 6),
             rep(NA, 10),
             0.098,
             -0.302,
             rep(NA, 6),
             rep(NA, 10),
             NA),
  FDR = c(0.732,
          rep(0.046, 6),
          rep(NA, 10),
          0.643,
          0.05,
          rep(NA, 6),
          rep(NA, 10),
          NA),
  groups = c(rep("ATP6V1E2", 18), rep("ADAP1", 18))
)

p1 <- dt_1 %>%
  group_by(groups) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = log2FC)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~groups) +
  theme_void() 

dt_1$is.sig <- case_when(
  abs(dt_1$log2FC) >= log2(1.1) & dt_1$FDR <= 0.05 ~ "Significant",
  is.na(dt_1$log2FC) ~ "Not expressed",
  TRUE ~ "Not significant"
)

dt_2 <- dt_1[dt_1$is.sig=="Significant",]
p2 <- dt_2 %>%
  group_by(groups) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = is.sig)) +
  facet_wrap(~groups) +
  theme_void()

pdf("brain_map_log2FC_ADAP1_ATP6V1E2.pdf", width = 6, height = 4)
grid.arrange(p1, p2, nrow = 2)
dev.off()

###############################################################################
setwd("/Users/phammaitam/Documents/research_projects/MD/Tinnitus-MCD/output/")

library(ggplot2)
library(ggseg)
library(gridExtra)

# Gene expression
dt_1 <- tibble(
  region = rep(c("entorhinal", 
                 "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
                 "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medical orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis", "precentral gyrus",
                 "parahippocampal",
                 "precentral", "paracentral",
                 "postcentral",
                 "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"), 2),
  log2FC = c(-0.061,
             rep(0.483, 6),
             rep(NA, 10),
             0.098,
             rep(NA,2),
             NA,
             rep(NA,4),
             -0.302,
             rep(NA, 6),
             rep(NA, 10),
             -0.268,
             rep(-0.145, 2),
             -0.176,
             rep(-0.205, 4)),
  FDR = c(0.732,
          rep(0.046, 6),
          rep(NA, 10),
          0.643,
          rep(NA,2),
          NA,
          rep(NA,4),
          0.011,
          rep(NA, 6),
          rep(NA, 10),
          0.008,
          rep(0.12, 2),
          0.008,
          rep(0.004, 4)),
  groups = c(rep("ATP6V1E2", 25), rep("ADAP1", 25))
)

p1 <- dt_1 %>%
  group_by(groups) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = log2FC)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~groups) +
  theme_void() 

dt_1$is.sig <- case_when(
  abs(dt_1$log2FC) >= log2(1.1) & dt_1$FDR <= 0.05 ~ "Significant",
  is.na(dt_1$log2FC) ~ "Not expressed",
  TRUE ~ "Not significant"
)

dt_2 <- dt_1[dt_1$is.sig=="Significant",]
p2 <- dt_2 %>%
  group_by(groups) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = is.sig)) +
  facet_wrap(~groups) +
  theme_void()

dt_3 <- data.frame(region=c("entorhinal", 
                            "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
                            "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medical orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis", "precentral gyrus",
                            "parahippocampal",
                            "precentral", "paracentral",
                            "postcentral",
                            "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"),
                   cat=c("Entorhinal Cortex (1,2)",
                         rep("Temporal Cortex (1)", 6),
                         rep("Frontal Cortex (1)", 10),
                         "Hippocampus (1,2)",
                         rep("Motor Cortex (2)",2),
                         "Sensory Cortex (2)",
                         rep("Cingulate Gyrus (2)", 4)))

p3 <- dt_3 %>%
 ggplot() +
 geom_brain(atlas = dk, 
            position = position_brain(hemi ~ side),
            aes(fill = cat)) +
 theme_void()
  
pdf("brain_map_log2FC_ADAP1_ATP6V1E2.pdf", width = 6, height = 6)
grid.arrange(
  grobs = list(p3, p1, p2),
  layout_matrix = rbind(
    c(1, 1),
    c(2,2),
    c(3,3)
  ),
)
dev.off()

################################################################################
dt_1 <- read.delim("brain_map.tsv", header=T, na.strings=c("."))

dt_all <- data.frame(region=c("entorhinal", 
                            "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
                            "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medial orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis",
                            "parahippocampal",
                            "precentral", "paracentral",
                            "postcentral",
                            "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"),
                   cat=c("Entorhinal Cortex",
                         rep("Temporal Cortex", 6),
                         rep("Frontal Cortex", 9),
                         "Hippocampus",
                         rep("Motor Cortex",2),
                         "Sensory Cortex",
                         rep("Cingulate Gyrus", 4)))


p_1 <- dt_all %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = cat)) +
  theme_void()

dt_2 <- dt_1 %>%
  group_by(Comparison, Brain.region) %>%
  summarise(ave.Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop")
dt_2$Comparison <- factor(dt_2$Comparison)
dt_2$Brain.region <- factor(dt_2$Brain.region)
dt_3 <- merge(dt_2, dt_all, by.y="cat", by.x="Brain.region", all.y=TRUE)

p_2 <- dt_3 %>%
  group_by(Comparison) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = ave.Log2FC)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~Comparison) +
  theme_void() 

dt_4 <- unique(dt_1[,c(1,2,3)])
dt_5 <- dt_4 %>%
  group_by(Comparison, Brain.region) %>%
  summarise(count=n(), .groups="drop")
dt_5 <- merge(dt_5, dt_all, by.y="cat", by.x="Brain.region", all.y=TRUE)

p_3 <- dt_5 %>%
  group_by(Comparison) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = count)) +
  scale_fill_gradient2(high = "cyan") +
  facet_wrap(~Comparison) +
  theme_void() 

matrix(x, nrow = 2, byrow=T)

#############################################################################
# prkra
dt_all <- data.frame(region=c("entorhinal", 
                              "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
                              "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medial orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis",
                              "parahippocampal",
                              "precentral", "paracentral",
                              "postcentral",
                              "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"),
                     cat=c("Entorhinal Cortex",
                           rep("Temporal Cortex", 6),
                           rep("Frontal Cortex", 9),
                           "Hippocampus",
                           rep("Motor Cortex",2),
                           "Sensory Cortex",
                           rep("Cingulate Gyrus", 4)))

# RNA-seq expression
dt_rna <- data.frame(Brain.region=c("Entorhinal Cortex", 
                              "Temporal Cortex",
                              "Frontal Cortex",
                              "Hippocampus",
                              "Motor Cortex",
                              "Sensory Cortex",
                              "Cingulate Gyrus"),
                     tpm=c(10.2,
                           7.5+7.1+6.1+6.1+9.5+4.4+8.4+6.4,
                           7.9+6.2+9.2+6.6+7.7+7+7.3+7.6,
                           10.2,
                           9.2+3.1+7.7+11.2+8.1,
                           6.2+9.7+6.4+5.1,
                           7+6.2+4.4+7.1+5.2+6.2))
dt_rna <- merge(dt_rna, dt_all, by.x="Brain.region", by.y="cat", all.y=TRUE)

p_1 <- dt_rna %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = tpm)) +
  scale_fill_gradient2(low = "#56B1F7", high = "#132B43") +
  theme_void()

# Protein expression change in AD vs control
dt <- read.delim("brain_map.tsv", header=T, na.strings=c("."))
dt_pr <- dt[dt$Gene=="PRKRA",]

dt_2 <- dt_pr %>%
  group_by(Brain.region) %>%
  summarise(ave.Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop")
dt_2$Brain.region <- factor(dt_2$Brain.region)
dt_3 <- merge(dt_2, dt_all, by.y="cat", by.x="Brain.region", all.y=TRUE)

p_2 <- dt_3 %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = ave.Log2FC)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_void()

pdf("brain_map_PRKRA.pdf", width = 4, height = 6)
grid.arrange(grobs = list(p_1, p_2), nrow=2, ncol=1)
dev.off()

#################################################################################
# eif4g1
dt_all <- data.frame(region=c("entorhinal", 
                              "inferior temporal", "superior temporal", "temporal pole", "transverse temporal", "middle temporal", "bankssts",
                              "caudal middle frontal", "frontal pole", "lateral orbitofrontal", "rostral middle frontal", "superior frontal", "medial orbitofrontal", "pars opercularis", "pars orbitalis", "pars triangularis",
                              "parahippocampal",
                              "precentral", "paracentral",
                              "postcentral",
                              "rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"),
                     cat=c("Entorhinal Cortex",
                           rep("Temporal Cortex", 6),
                           rep("Frontal Cortex", 9),
                           "Hippocampus",
                           rep("Motor Cortex",2),
                           "Sensory Cortex",
                           rep("Cingulate Gyrus", 4)))

# RNA-seq expression
dt_rna <- data.frame(Brain.region=c("Entorhinal Cortex", 
                                    "Temporal Cortex",
                                    "Frontal Cortex",
                                    "Hippocampus",
                                    "Motor Cortex",
                                    "Sensory Cortex",
                                    "Cingulate Gyrus"),
                     tpm=c(118.1,
                           98.3+90.1+101.3+97.9+105.9+93.1+90.7+91.7,
                           90.1+87.1+99.5+95.4+97.6+129+129.5+105.3,
                           118.1,
                           97.9+110.9+84+84.6+93.8,
                           105+93+94.5+86.1,
                           96.9+98.8+100.2+100+103.3+92.5))
dt_rna <- merge(dt_rna, dt_all, by.x="Brain.region", by.y="cat", all.y=TRUE)

p_1 <- dt_rna %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = tpm)) +
  scale_fill_gradient2(low = "#56B1F7", high = "#132B43") +
  theme_void()

# Protein expression change in AD vs control
dt <- read.delim("brain_map.tsv", header=T, na.strings=c("."))
dt_pr <- dt[dt$Gene=="EIF4G1",]

dt_2 <- dt_pr %>%
  group_by(Brain.region) %>%
  summarise(ave.Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop")
dt_2$Brain.region <- factor(dt_2$Brain.region)
dt_3 <- merge(dt_2, dt_all, by.y="cat", by.x="Brain.region", all.y=TRUE)

p_2 <- dt_3 %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = ave.Log2FC)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_void()

pdf("brain_map_EIF4G1.pdf", width = 4, height = 6)
grid.arrange(grobs = list(p_1, p_2), nrow=2, ncol=1)
dev.off()