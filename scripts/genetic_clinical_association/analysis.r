# revised 17/10/2025

#################################################################################
### Clinical scores and APOE risk analysis
#################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(cowplot)
PRJ_DIR <- "/home/mdnl/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1"

GT <- read.table(paste(PRJ_DIR, "APOE_risk/APOE_GT_rs429358_rs7412_LF.tsv", sep="/"), header=T, na.strings = c(".", "NA", "N/A"))
GT$rs429358 <- gsub("\\|", "/", GT$rs429358)
GT$rs7412 <-gsub("\\|", "/", GT$rs7412)

GT <- GT %>%
  mutate(
    APOE_risk = case_when(
      rs429358 == "0/0" & rs7412 == "1/1" ~ "ε2/ε2",
      rs429358 == "0/0" & rs7412 == "0/1" ~ "ε2/ε3",
      rs429358 == "0/1" & rs7412 == "0/1" ~ "ε2/ε4",
      rs429358 == "0/0" & rs7412 == "0/0" ~ "ε3/ε3",
      rs429358 == "0/1" & rs7412 == "0/0" ~ "ε3/ε4",
      rs429358 == "1/1" & rs7412 == "0/0" ~ "ε4/ε4",
      TRUE ~ NA_character_
    )
  )
GT$APOE_risk <- as.factor(GT$APOE_risk)

THI_GUEF <- read.table(paste(PRJ_DIR, "APOE_risk/THI_GUEF_score.tsv", sep="/"), header=T, na.strings = c(".", "NA", "N/A"))

# Select THI scores
THI_sub <- THI_GUEF[THI_GUEF$visit_type %in% c("baseline", "screening"), c(1,2,5,8:10)]

# Convert visit_day to Date if not already
THI_sub$visit_day <- as.Date(THI_sub$visit_day)

# Sort data: by patient_code, then descending visit_day
THI_sub <- THI_sub[order(THI_sub$SAMPLE_ID, THI_sub$visit_type, -as.numeric(THI_sub$visit_day)), ]

# Split by patient_code
patients <- split(THI_sub, THI_sub$SAMPLE_ID)

# Function to apply the logic for one patient
get_best_thi <- function(pdata) {
  # Get all baselines and screenings
  baseline_rows <- pdata[pdata$visit_type == "baseline", ]
  screening_rows <- pdata[pdata$visit_type == "screening", ]

  n_baseline <- nrow(baseline_rows)
  
  if (n_baseline == 1) {
      if (!is.na(baseline_rows$thi_score)) {
          return(baseline_rows)
      } else {
      # return the most recent screening with non-NA thi
      non_na_screening <- screening_rows[!is.na(screening_rows$thi_score), ]
          if (nrow(non_na_screening) > 0) {
              return(non_na_screening[1,])  # already sorted descending
          }
      }
  }

  if (n_baseline >= 2) {
      non_na_baseline <- baseline_rows[!is.na(baseline_rows$thi_score), ]
      if (nrow(non_na_baseline) == 1) {
          return(non_na_baseline)
      } else if (nrow(non_na_baseline) >= 2) {
          return(non_na_baseline[1,])  # most recent
      }
  }

  # fallback: return most recent screening with non-NA THI
  non_na_screening <- screening_rows[!is.na(screening_rows$thi_score), ]
  if (nrow(non_na_screening) > 0) {
      return(non_na_screening[1,])
  }

  return(pdata[1,])
}

best_thi_scores <- sapply(patients, get_best_thi)

# Convert to a data frame
best_thi <- data.frame(t(best_thi_scores))                    
best_thi <- best_thi[,-c(4:5)]

# GUEF score
GUEF_sub <- THI_GUEF[THI_GUEF$visit_type %in% c("baseline", "screening"), c(1,2,5,8,9,11)]

# Convert visit_day to Date if not already
GUEF_sub$visit_day <- as.Date(GUEF_sub$visit_day)

# Sort data: by patient_code, then descending visit_day
GUEF_sub <- GUEF_sub[order(GUEF_sub$SAMPLE_ID, GUEF_sub$visit_type, -as.numeric(GUEF_sub$visit_day)), ]

# Split by patient_code
patients <- split(GUEF_sub, GUEF_sub$SAMPLE_ID)

# Function to apply the logic for one patient
get_best_guef <- function(pdata) {
  # Get all baselines and screenings
  baseline_rows <- pdata[pdata$visit_type == "baseline", ]
  screening_rows <- pdata[pdata$visit_type == "screening", ]

  n_baseline <- nrow(baseline_rows)
  
  if (n_baseline == 1) {
      if (!is.na(baseline_rows$guef_score)) {
          return(baseline_rows)
      } else {
      # return the most recent screening with non-NA thi
      non_na_screening <- screening_rows[!is.na(screening_rows$guef_score), ]
          if (nrow(non_na_screening) > 0) {
              return(non_na_screening[1,])  # already sorted descending
          }
      }
  }

  if (n_baseline >= 2) {
      non_na_baseline <- baseline_rows[!is.na(baseline_rows$guef_score), ]
      if (nrow(non_na_baseline) == 1) {
          return(non_na_baseline)
      } else if (nrow(non_na_baseline) >= 2) {
          return(non_na_baseline[1,])  # most recent
      }
  }

  # fallback: return most recent screening with non-NA THI
  non_na_screening <- screening_rows[!is.na(screening_rows$guef_score), ]
  if (nrow(non_na_screening) > 0) {
      return(non_na_screening[1,])
  }

  return(pdata[1,])
}

best_guef_scores <- sapply(patients, get_best_guef)

# Convert to a data frame
best_guef <- data.frame(t(best_guef_scores))                    
best_guef <- best_guef[,c(2,6)]          

# MOCA score
MOCA <- read.csv(paste(PRJ_DIR, "APOE_risk/MOCA_score.csv", sep="/"), header=T)
MOCA <- MOCA %>%                                                                        
  group_by(sample_ID) %>%
  summarise(MoCA.Score = mean(MoCA.Score, na.rm = TRUE))

# HL score
HL <- read.csv(paste(PRJ_DIR, "APOE_risk/HL_UNITI.csv", sep="/"), header=T, na.strings = c(".", "NA", "N/A"))
HL <- HL[, c(1,8,14:16,17:26)]
colnames(HL) <- c("SAMPLE_ID", "HL_sex", "HL_thi_score", "HL_guef_score", "HL_PHQ9_score",
                  "HL_500Hz_right", 
                  "HL_1kHz_right",
                  "HL_2kHz_right",
                  "HL_500Hz_left",
                  "HL_1kHz_left",
                  "HL_2kHz_left",
                  "HL_4kHz_right",
                  "HL_8kHz_right",
                  "HL_4kHz_left",
                  "HL_8kHz_left")
HL[6:15][HL[6:15] > 120] <- NA

APOE_THI <- merge(GT, best_thi, by = "SAMPLE_ID")
APOE_THI_GUEF <- merge(APOE_THI, best_guef, by = "SAMPLE_ID")
APOE_THI_GUEF_MOCA <- merge(APOE_THI_GUEF, MOCA, by.x = "SAMPLE_ID", by.y="sample_ID", all.x=TRUE)
APOE_THI_GUEF_MOCA <- merge(APOE_THI_GUEF_MOCA, HL, by = "SAMPLE_ID", all.x=TRUE)

# # Manually create samples of THI 5 grades for GBA analysis
# dir.create("GBA_THI", showWarnings = FALSE, recursive = TRUE)
# thi_g1_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score<=16, "SAMPLE_ID"]
# write.table(as.data.frame(thi_g1_samples), file="GBA_THI/thi_g1_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_g2_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score<=36 & APOE_THI_GUEF_MOCA$thi_score >=18, "SAMPLE_ID"]
# write.table(as.data.frame(thi_g2_samples), file="GBA_THI/thi_g2_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_g3_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score<=56 & APOE_THI_GUEF_MOCA$thi_score >=38, "SAMPLE_ID"]
# write.table(as.data.frame(thi_g3_samples), file="GBA_THI/thi_g3_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_g4_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score<=76 & APOE_THI_GUEF_MOCA$thi_score >=58, "SAMPLE_ID"]
# write.table(as.data.frame(thi_g4_samples), file="GBA_THI/thi_g4_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_g5_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score >=78, "SAMPLE_ID"]
# write.table(as.data.frame(thi_g5_samples), file="GBA_THI/thi_g5_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_mild_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score<=36, "SAMPLE_ID"]
# write.table(as.data.frame(thi_mild_samples), file="GBA_THI/thi_mild_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_severe_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$thi_score) & APOE_THI_GUEF_MOCA$thi_score>=58, "SAMPLE_ID"]
# write.table(as.data.frame(thi_severe_samples), file="GBA_THI/thi_severe_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# # Manually create samples of THI 5 grades for GBA analysis
# dir.create("GBA_THI", showWarnings = FALSE, recursive = TRUE)
# guef_g1_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$guef_score) & APOE_THI_GUEF_MOCA$guef_score<=10, "SAMPLE_ID"]
# write.table(as.data.frame(guef_g1_samples), file="GBA_GUEF/guef_g1_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# guef_g2_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$guef_score) & APOE_THI_GUEF_MOCA$guef_score<=17 & APOE_THI_GUEF_MOCA$guef_score >=11, "SAMPLE_ID"]
# write.table(as.data.frame(guef_g2_samples), file="GBA_GUEF/guef_g2_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# guef_g3_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$guef_score) & APOE_THI_GUEF_MOCA$guef_score<=25 & APOE_THI_GUEF_MOCA$guef_score >=18, "SAMPLE_ID"]
# write.table(as.data.frame(guef_g3_samples), file="GBA_GUEF/guef_g3_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# guef_g4_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$guef_score) & APOE_THI_GUEF_MOCA$guef_score<=45 & APOE_THI_GUEF_MOCA$guef_score >=26, "SAMPLE_ID"]
# write.table(as.data.frame(guef_g4_samples), file="GBA_GUEF/guef_g4_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# guef_severe_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$guef_score) & APOE_THI_GUEF_MOCA$guef_score >=18, "SAMPLE_ID"]
# write.table(as.data.frame(guef_severe_samples), file="GBA_GUEF/guef_severe_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# # Manually create samples of THI 5 grades for GBA analysis
# dir.create("GBA_PHQ9", showWarnings = FALSE, recursive = TRUE)
# phq9_g1_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$HL_PHQ9_score) & APOE_THI_GUEF_MOCA$HL_PHQ9_score<=4, "SAMPLE_ID"]
# write.table(as.data.frame(phq9_g1_samples), file="GBA_PHQ9/phq9_g1_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# phq9_g2_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$HL_PHQ9_score) & APOE_THI_GUEF_MOCA$HL_PHQ9_score<=9 & APOE_THI_GUEF_MOCA$HL_PHQ9_score>=5, "SAMPLE_ID"]
# write.table(as.data.frame(phq9_g2_samples), file="GBA_PHQ9/phq9_g2_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# phq9_g3_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$HL_PHQ9_score) & APOE_THI_GUEF_MOCA$HL_PHQ9_score<=14 & APOE_THI_GUEF_MOCA$HL_PHQ9_score>=10, "SAMPLE_ID"]
# write.table(as.data.frame(phq9_g3_samples), file="GBA_PHQ9/phq9_g3_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# phq9_g4_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$HL_PHQ9_score) & APOE_THI_GUEF_MOCA$HL_PHQ9_score<=27 & APOE_THI_GUEF_MOCA$HL_PHQ9_score>=15, "SAMPLE_ID"]
# write.table(as.data.frame(phq9_g4_samples), file="GBA_PHQ9/phq9_g4_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# phq9_severe_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$HL_PHQ9_score) & APOE_THI_GUEF_MOCA$HL_PHQ9_score>=10, "SAMPLE_ID"]
# write.table(as.data.frame(phq9_severe_samples), file="GBA_PHQ9/phq9_severe_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# # Manually create samples of MOCA 2 grades for GBA analysis
# dir.create("GBA_revised", showWarnings = FALSE, recursive = TRUE)
# moca_normal_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.Score) & APOE_THI_GUEF_MOCA$MoCA.Score>25, "SAMPLE_ID"]
# write.table(as.data.frame(moca_normal_samples), file="GBA_revised/moca_normal_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Start analyze APOE_THI_GUEF_MOCA
APOE_THI_GUEF_MOCA$APOE_risk <- as.character(APOE_THI_GUEF_MOCA$APOE_risk)
APOE_THI_GUEF_MOCA$thi_score <- as.numeric(APOE_THI_GUEF_MOCA$thi_score)
APOE_THI_GUEF_MOCA$guef_score <- as.numeric(APOE_THI_GUEF_MOCA$guef_score)
APOE_THI_GUEF_MOCA$MoCA.Score <- as.numeric(APOE_THI_GUEF_MOCA$MoCA.Score)
APOE_THI_GUEF_MOCA$age <- as.numeric(APOE_THI_GUEF_MOCA$age)

APOE_THI_GUEF_MOCA$APOE_risk_1 <- APOE_THI_GUEF_MOCA$APOE_risk
APOE_THI_GUEF_MOCA$APOE_risk_1[APOE_THI_GUEF_MOCA$APOE_risk %in% c("ε2/ε3", "ε3/ε3")] <- "Non-APOE E4"
APOE_THI_GUEF_MOCA$APOE_risk_1[APOE_THI_GUEF_MOCA$APOE_risk %in% c("ε2/ε4", "ε3/ε4", "ε4/ε4")] <- "APOE E4"
APOE_THI_GUEF_MOCA$APOE_risk_1 <- factor(APOE_THI_GUEF_MOCA$APOE_risk_1, levels = c("Non-APOE E4", "APOE E4"))

thi_q4 <- quantile(APOE_THI_GUEF_MOCA$thi_score, probs=0.75, na.rm=T)[[1]]
thi_q1 <- quantile(APOE_THI_GUEF_MOCA$thi_score, probs=0.25, na.rm=T)[[1]]
APOE_THI_GUEF_MOCA$tinnitus.level <- ifelse(APOE_THI_GUEF_MOCA$thi_score >= thi_q4, "Q4",
                                  ifelse(APOE_THI_GUEF_MOCA$thi_score <= thi_q1, "Q1", "normal"))
APOE_THI_GUEF_MOCA$tinnitus.level <- factor(APOE_THI_GUEF_MOCA$tinnitus.level, levels = c("Q1", "normal", "Q4"))

# Manually create samples of THI_Q4 for GBA analysis
# dir.create("GBA_THI", showWarnings = FALSE, recursive = TRUE)
# thi_q4_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$tinnitus.level) & APOE_THI_GUEF_MOCA$tinnitus.level == "Q4", "SAMPLE_ID"]
# write.table(as.data.frame(thi_q4_samples), file="GBA_THI/thi_q4_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# thi_q1_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$tinnitus.level) & APOE_THI_GUEF_MOCA$tinnitus.level == "Q1", "SAMPLE_ID"]
# write.table(as.data.frame(thi_q1_samples), file="GBA_THI/thi_q1_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

guef_q4 <- quantile(APOE_THI_GUEF_MOCA$guef_score, probs=0.75, na.rm=T)[[1]]
guef_q1 <- quantile(APOE_THI_GUEF_MOCA$guef_score, probs=0.25, na.rm=T)[[1]]
APOE_THI_GUEF_MOCA$hyperacusis.level <- ifelse(APOE_THI_GUEF_MOCA$guef_score >= guef_q4, "Q4",
                                  ifelse(APOE_THI_GUEF_MOCA$guef_score <= guef_q1, "Q1", "normal"))
APOE_THI_GUEF_MOCA$hyperacusis.level <- factor(APOE_THI_GUEF_MOCA$hyperacusis.level, levels = c("Q1", "normal", "Q4"))

# # Manually create samples of GUEF_Q4 for GBA analysis
# dir.create("GBA_GUEF", showWarnings = FALSE, recursive = TRUE)
# guef_q4_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$hyperacusis.level) & APOE_THI_GUEF_MOCA$hyperacusis.level == "Q4", "SAMPLE_ID"]
# write.table(as.data.frame(guef_q4_samples), file="GBA_GUEF/guef_q4_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# guef_q1_samples <- APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$hyperacusis.level) & APOE_THI_GUEF_MOCA$hyperacusis.level == "Q1", "SAMPLE_ID"]
# write.table(as.data.frame(guef_q1_samples), file="GBA_GUEF/guef_q1_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

MoCA_q1 <- 25
APOE_THI_GUEF_MOCA$MoCA.level <- ifelse(APOE_THI_GUEF_MOCA$MoCA.Score <= MoCA_q1, "Q1", "normal")
APOE_THI_GUEF_MOCA$MoCA.level <- factor(APOE_THI_GUEF_MOCA$MoCA.level, levels = c("Q1", "normal"))

APOE_THI_GUEF_MOCA$age.level <- ifelse(APOE_THI_GUEF_MOCA$age < 60, "age < 60", "age >= 60")
APOE_THI_GUEF_MOCA$age.level <- factor(APOE_THI_GUEF_MOCA$age.level, levels = c("age < 60", "age >= 60"))

APOE_THI_GUEF_MOCA$PTA_right <- rowMeans(APOE_THI_GUEF_MOCA[, c("HL_500Hz_right", "HL_1kHz_right", "HL_2kHz_right")], na.rm=TRUE)
APOE_THI_GUEF_MOCA$PTA_right.HL <- ifelse(APOE_THI_GUEF_MOCA$PTA_right > 25, "HL", "non-HL")
APOE_THI_GUEF_MOCA$PTA_left <- rowMeans(APOE_THI_GUEF_MOCA[, c("HL_500Hz_left", "HL_1kHz_left", "HL_2kHz_left")], na.rm=TRUE)
APOE_THI_GUEF_MOCA$PTA_left.HL <- ifelse(APOE_THI_GUEF_MOCA$PTA_left > 25, "HL", "non-HL")

APOE_THI_GUEF_MOCA$PTA_bilateral.HL <- ifelse(APOE_THI_GUEF_MOCA$PTA_right.HL == "HL" & APOE_THI_GUEF_MOCA$PTA_left.HL == "HL", "bilateral-HL",
                                        ifelse(APOE_THI_GUEF_MOCA$PTA_right.HL == "HL" | APOE_THI_GUEF_MOCA$PTA_left.HL == "HL", "unilateral-HL", "non-HL"))
APOE_THI_GUEF_MOCA$PTA_bilateral.HL <- factor(APOE_THI_GUEF_MOCA$PTA_bilateral.HL, levels = c("non-HL", "unilateral-HL", "bilateral-HL"))

APOE_THI_GUEF_MOCA$HF_right.HL <- ifelse(APOE_THI_GUEF_MOCA$HL_8kHz_right > 25 & (APOE_THI_GUEF_MOCA$HL_8kHz_right - APOE_THI_GUEF_MOCA$PTA_right) > 15, "HL", "non-HL")
APOE_THI_GUEF_MOCA$HF_left.HL <- ifelse(APOE_THI_GUEF_MOCA$HL_8kHz_left > 25 & (APOE_THI_GUEF_MOCA$HL_8kHz_left - APOE_THI_GUEF_MOCA$PTA_left) > 15, "HL", "non-HL")
APOE_THI_GUEF_MOCA$HF_bilateral.HL <- ifelse(APOE_THI_GUEF_MOCA$HF_right.HL == "HL" & APOE_THI_GUEF_MOCA$HF_left.HL == "HL", "bilateral-HL",
                                        ifelse(APOE_THI_GUEF_MOCA$HF_right.HL == "HL" | APOE_THI_GUEF_MOCA$HF_left.HL == "HL", "unilateral-HL", "non-HL"))
APOE_THI_GUEF_MOCA$HF_bilateral.HL <- factor(APOE_THI_GUEF_MOCA$HF_bilateral.HL, levels = c("non-HL", "unilateral-HL", "bilateral-HL"))
APOE_THI_GUEF_MOCA$HF.HL <- rowMeans(APOE_THI_GUEF_MOCA[, c("HL_4kHz_right", "HL_8kHz_right", "HL_4kHz_left", "HL_8kHz_left")], na.rm=TRUE)
APOE_THI_GUEF_MOCA$PTA <- rowMeans(APOE_THI_GUEF_MOCA[, c("PTA_right", "PTA_left")], na.rm=TRUE) 

# MOCA distribution
moca_q1 <- quantile(APOE_THI_GUEF_MOCA$MoCA.Score, 0.25, na.rm = TRUE)
moca_q4 <- quantile(APOE_THI_GUEF_MOCA$MoCA.Score, 0.75, na.rm = TRUE)

# Histogram for age vs MCI/non-MCI groups
png(paste(PRJ_DIR, "figures/MCI_vs_nonMCI_age_dist.png", sep="/"), width=5*300, height=4*300, res=300)
ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level), ], aes(x = age, fill = MoCA.level)) +
  geom_density(alpha = 0.4) +
  labs(x = "Age", 
    y = "Density", 
    title = "Distribution of inclusion age in MCI vs non-MCI groups") +
  scale_fill_discrete(
    name = "Cognitive status",
    labels = c(normal = "Non-MCI", Q1 = "MCI")
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, hjust = 0.5, color = "black")
  )
dev.off()

# Draw histogram
sites <- substr(APOE_THI_GUEF_MOCA$SAMPLE_ID, 1, 2)
df <- as.data.frame(table(sites))
colnames(df) <- c("Site", "Count")
p1 <- ggplot(df, aes(x = Site, y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample distribution across clinical centers",
       x = "Clinical center",
       y = "Number of samples") +
  theme_minimal()

p2 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = MoCA.Score)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "gray85", color = "black") +
  geom_density(color = "darkblue", size = 1.2, alpha = 0.7) +
  geom_vline(xintercept = moca_q1, color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = moca_q4, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = moca_q1-0.5, y = 0.17, label = "Q1", color = "blue", vjust = -1, fontface = "bold", size=4) +
  annotate("text", x = moca_q4-0.5, y = 0.17, label = "Q4", color = "red", vjust = -1, fontface = "bold", size=4) +
  labs(x = "MoCA Score",
      title = "Distribution of MoCA scores (n = 276)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 12, hjust = 0.5, color = "black")
  )

png(paste(PRJ_DIR, "figures/MOCA_score_distribution.png", sep="/"), width = 8*300, height = 4*300, res=300)
plot_grid(p1, p2, nrow=1, labels = c("A", "B"))
dev.off()

# Statistical testing
# mod <- summary(lm(MoCA.Score ~ thi_score, data = APOE_THI_GUEF_MOCA))
# data.frame(
#   r = sign(coef(mod)[2]) * sqrt(out$r.squared),
#   R2 = out$r.squared,
#   p_value = out$coefficients[2, 4]
# )
library(Hmisc)
library(reshape2)

df <- list()
df[["all"]] <- APOE_THI_GUEF_MOCA
df[["MCI"]] <- APOE_THI_GUEF_MOCA[APOE_THI_GUEF_MOCA$MoCA.Score <= 25, ]
df[["Non-MCI"]] <- APOE_THI_GUEF_MOCA[APOE_THI_GUEF_MOCA$MoCA.Score > 25, ]
for (group in names(df)) {

  # correlation matrix
  print(group)
  cor.mat <- rcorr(as.matrix(df[[group]][, c("MoCA.Score", "thi_score", "guef_score", "PTA", "HF.HL", "HL_PHQ9_score")]), type = "pearson")
  res_df <- melt(cor.mat$r, varnames = c("Var1", "Var2"), value.name = "cor")
  p_df <- melt(cor.mat$P, varnames = c("Var1", "Var2"), value.name = "pval")

  # Merge correlation and p-value data
  res_df <- merge(res_df, p_df, by = c("Var1", "Var2"))
  res_df <- subset(res_df, Var1 != Var2)
  res_df <- res_df[!duplicated(t(apply(res_df[,1:2], 1, sort))), ]
  res_df$pval_adj <- p.adjust(res_df$pval, method = "BH")
  print(res_df)
}

# Significant gene burden
LoF_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.LoF.GBA.v3.csv", header=T)
missense_dt <- read.csv("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.missense.GBA.v3.csv", header=T)
synap_gba <- read.table("GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.GBA_genes.v3.txt", header=TRUE)
synap_genes <- synap_gba$GBA_gene
BD_genes <- c("PRKRA", "GRAMD1A", "ENDOG", "PIK3C2G", "IGFN1", "EIF4G1", "SCN5A", "ERBB3", "TBCC", "PRX", "ANKRD52", "ANO2")
AD_genes <- c("PRKRA", "EIF4G1")
tinnitus_genes <- c("GRAMD1A", "ANKRD52", "ENDOG", "NCKAP5", "LGALS12", "PSKH1", "KCNK12", "UPK2", "CDHR4", "PLEKHA3")
tinnitus_synap_genes <- c("GRAMD1A", "ANKRD52", "ENDOG")
ns_tinnitus_genes <- c("PRKRA", "IGFN1", "PIK3C2G", "PRX", "HK3", "ATP6V1E2", "SLC12A8", "DLX4", "ENTPD4", "TTC6", "CAPN12", "ERMARD", "OR5AU1")
ns_tinnitus_synap_genes <- c("PRKRA", "IGFN1", "PIK3C2G", "PRX", "HK3", "ATP6V1E2")
hyperacusis_genes <- c("ANKRD52", "SCN5A", "NCKAP5", "RBM23")
hyperacusis_synap_genes <- c("ANKRD52", "SCN5A")

LoF_samples <- LoF_dt[, c("SYMBOL", "samples_het", "samples_hom")]
missense_samples <- missense_dt[, c("SYMBOL", "samples_het", "samples_hom")]
sample_data <- rbind(LoF_samples, missense_samples)

# SV
SV_001 <- read.csv("SV_CNV_output/SV_GOI_overlap_revised.TIDDIT.001.csv", header=T, na.strings = c(".", "NA", "N/A"))
SV_011 <- read.csv("SV_CNV_output/SV_GOI_overlap_revised.v3.csv", header=T, na.strings = c(".", "NA", "N/A"))

target_SV <- data.frame(gene = c("PRKRA", "TNKS1BP1", "ATP6V1E2", "PIK3C2G"), 
                        SV_pos = c("2_178444502_178447509_DEL", "11_57304632_57304828_DEL", 
                                   "2_46524048_46524172_DEL", "12_18570336_18572160_DEL"),
                        samples_001 = rep(NA,4),
                        samples_011 = rep(NA,4))


for (g in target_SV$gene) {
  sv_pos <- target_SV$SV_pos[target_SV$gene == g]

  idx <- which(SV_001$Gene == g)
  samples_001 <- colnames(SV_001)[grep(sv_pos, SV_001[idx, ])]
  if (length(samples_001) == 0) {
    samples_001 <- NA
  } else {
    samples_001 <- paste(samples_001, collapse = "|")
  }
  target_SV$samples_001[target_SV$gene == g] <- samples_001

  idx <- which(SV_011$Gene == g)
  samples_011 <- colnames(SV_011)[grep(sv_pos, SV_011[idx, ])]
  if (length(samples_011) == 0) {
    samples_011 <- NA
  } else {
    samples_011 <- paste(samples_011, collapse = "|")
  }
  target_SV$samples_011[target_SV$gene == g] <- samples_011

}


total_samples <- list()

for (row in 1:nrow(sample_data)) {
    print(sample_data$SYMBOL[row])
    s1 <- unique(unlist(strsplit(sample_data$samples_het[row], "[|;]")))
    s2 <- unique(unlist(strsplit(sample_data$samples_hom[row], "[|;]")))
    s1 <- s1[!is.na(s1) & s1 != "NA"]
    s2 <- s2[!is.na(s2) & s2 != "NA"]
    if (sample_data$SYMBOL[row] %in% BD_genes) {
      # dt_1 <- APOE_THI_GUEF_MOCA[APOE_THI_GUEF_MOCA$SAMPLE_ID %in% c(s1,s2), ]
      total_samples[["BD_genes"]] <- unique(c(total_samples[["BD_genes"]], s1, s2))
    }

    if (sample_data$SYMBOL[row] %in% synap_genes) {
      total_samples[["synap_genes"]] <- unique(c(total_samples[["synap_genes"]], s1, s2))
    }

    if (sample_data$SYMBOL[row] %in% AD_genes) {
      total_samples[["AD_genes"]] <- unique(c(total_samples[["AD_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% tinnitus_genes) {
      total_samples[["tinnitus_genes"]] <- unique(c(total_samples[["tinnitus_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% tinnitus_synap_genes) {
      total_samples[["tinnitus_synap_genes"]] <- unique(c(total_samples[["tinnitus_synap_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% ns_tinnitus_genes) {
      total_samples[["ns_tinnitus_genes"]] <- unique(c(total_samples[["ns_tinnitus_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% ns_tinnitus_synap_genes) {
      total_samples[["ns_tinnitus_synap_genes"]] <- unique(c(total_samples[["ns_tinnitus_synap_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% hyperacusis_genes) {
      total_samples[["hyperacusis_genes"]] <- unique(c(total_samples[["hyperacusis_genes"]], s1, s2))
    } 

    if (sample_data$SYMBOL[row] %in% hyperacusis_synap_genes) {
      total_samples[["hyperacusis_synap_genes"]] <- unique(c(total_samples[["hyperacusis_synap_genes"]], s1, s2))
    } 

    total_samples[["gba_genes"]] <- unique(c(total_samples[["gba_genes"]], s1, s2))
}


SV_genes <- target_SV$gene
for (row in 1:nrow(sample_data)) {
    print(sample_data$SYMBOL[row])
    if (sample_data$SYMBOL[row] %in% SV_genes) {
      s1 <- unique(unlist(strsplit(sample_data$samples_het[row], "[|;]")))
      s2 <- unique(unlist(strsplit(sample_data$samples_hom[row], "[|;]")))
      s3 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_001[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
      s4 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_011[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
      s1 <- s1[!is.na(s1) & s1 != "NA"]
      s2 <- s2[!is.na(s2) & s2 != "NA"]
      s3 <- s3[!is.na(s3) & s3 != "NA"]
      s4 <- s4[!is.na(s4) & s4 != "NA"]
      total_samples[["SV_genes"]] <- unique(c(total_samples[["SV_genes"]], s1, s2, s3, s4))
    }

    if (sample_data$SYMBOL[row] == "PRKRA") {
      s1 <- unique(unlist(strsplit(sample_data$samples_het[row], "[|;]")))
      s2 <- unique(unlist(strsplit(sample_data$samples_hom[row], "[|;]")))
      s3 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_001[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
      s4 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_011[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
      s1 <- s1[!is.na(s1) & s1 != "NA"]
      s2 <- s2[!is.na(s2) & s2 != "NA"]
      s3 <- s3[!is.na(s3) & s3 != "NA"]
      s4 <- s4[!is.na(s4) & s4 != "NA"]
      total_samples[["PRKRA"]] <- unique(c(total_samples[["PRKRA"]], s1, s2, s3, s4))
    }

}

# total_samples <- list()
# for (row in 1:nrow(sample_data)) {
#     print(sample_data$SYMBOL[row])
#     s1 <- unique(unlist(strsplit(sample_data$samples_het[row], "[|;]")))
#     s2 <- unique(unlist(strsplit(sample_data$samples_hom[row], "[|;]")))
#     if (sample_data$SYMBOL[row] %in% target_SV$gene) {
#       s3 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_001[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
#       s4 <- gsub("Sample", "", unlist(strsplit(target_SV$samples_011[target_SV$gene == sample_data$SYMBOL[row]], "\\|")))
#     }
#     s1 <- s1[!is.na(s1) & s1 != "NA"]
#     s2 <- s2[!is.na(s2) & s2 != "NA"]
#     s3 <- s3[!is.na(s3) & s3 != "NA"]
#     s4 <- s4[!is.na(s4) & s4 != "NA"]
#     if (sample_data$SYMBOL[row] %in% BD_genes) {
#       # dt_1 <- APOE_THI_GUEF_MOCA[APOE_THI_GUEF_MOCA$SAMPLE_ID %in% c(s1,s2), ]
#       total_samples[["BD_genes"]] <- unique(c(total_samples[["BD_genes"]], s1, s2, s3, s4))
#     }

#     if (sample_data$SYMBOL[row] %in% synap_genes) {
#       total_samples[["synap_genes"]] <- unique(c(total_samples[["synap_genes"]], s1, s2, s3, s4))
#     }

#     if (sample_data$SYMBOL[row] %in% AD_genes) {
#       total_samples[["AD_genes"]] <- unique(c(total_samples[["AD_genes"]], s1, s2, s3, s4))
#     } 

#     if (sample_data$SYMBOL[row] %in% tinnitus_genes) {
#       total_samples[["tinnitus_genes"]] <- unique(c(total_samples[["tinnitus_genes"]], s1, s2, s3, s4))
#     } 

#     total_samples[["gba_genes"]] <- unique(c(total_samples[["gba_genes"]], s1, s2, s3, s4))
# }

dt <- lapply(total_samples, function(x) APOE_THI_GUEF_MOCA[APOE_THI_GUEF_MOCA$SAMPLE_ID %in% x, ])
dt$PRKRA <- dt$PRKRA[dt$PRKRA$thi_score <= 56, ]

gene_types <- c("Brain disease genes", "Synaptic genes", "Alzheimer's disease genes", 
            "Severe tinnitus genes", "Severe tinnitus & synaptic genes","Non-severe tinnitus genes", "Non-severe tinnitus & synaptic genes", "All GBA genes", "SV genes", "PRKRA SNV/SV", "Severe hyperacusis genes", "Severe hyperacusis & synaptic genes")
names(gene_types) <- c("BD_genes", "synap_genes", "AD_genes", "tinnitus_genes", "tinnitus_synap_genes", "ns_tinnitus_genes", "ns_tinnitus_synap_genes", "gba_genes", "SV_genes", "PRKRA", "hyperacusis_genes", "hyperacusis_synap_genes")

y_max <- max(APOE_THI_GUEF_MOCA$thi_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$thi_score, na.rm=TRUE) - 10
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = MoCA.Score, y = thi_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "MoCA Score",
    y = "THI Score",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$guef_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$guef_score, na.rm=TRUE) - 10
p0.2 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = MoCA.Score, y = guef_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "MoCA Score",
    y = "GUF Score"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = MoCA.Score, y = thi_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +   
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "MoCA Score",
       y = "THI Score",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = MoCA.Score, y = guef_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  
   stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "MoCA Score",
       y = "GUF Score") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p2) <- names(dt)


# pdf("figures/THI_GUEF_MOCA_scatter_plot_v3.pdf", width=25, height=15)
png("figures/THI_GUEF_MOCA_scatter_plot_v3.png", width=25*300, height=15*300, res=300)
plot_grid(p0.1, p0.2, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]],
                        nrow=5, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", ""))
dev.off()


# scatter plot for MoCA score by HL level
y_max <- max(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) - 10
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = MoCA.Score, y = PTA, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "MoCA Score",
    y = "PTA (dB HL)",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) - 10
p0.2 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = MoCA.Score, y = HF.HL, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "MoCA Score",
    y = "High requency (dB HL)"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = MoCA.Score, y = PTA)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "MoCA Score",
       y = "PTA (dB HL)",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = MoCA.Score, y = HF.HL)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +
   stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "MoCA Score",
       y = "High frequency (dB HL)") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p2) <- names(dt)

# pdf("figures/HL_MOCA_scatter_plot_v2.pdf", width=25, height=15)
png("figures/HL_MOCA_scatter_plot_v2.png", width=25*300, height=15*300, res=300)
plot_grid(p0.1, p0.2, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]],
                        nrow=5, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", ""))
dev.off()

# scatter plot for MoCA score by depression level
y_max <- max(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) - 8
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = MoCA.Score, y = HL_PHQ9_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "MoCA Score",
    y = "PHQ-9 Score",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) - 8
p0.2 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = thi_score, y = HL_PHQ9_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "THI Score",
    y = "PHQ-9 Score"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) - 8
p0.3 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = guef_score, y = HL_PHQ9_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "GUF Score",
    y = "PHQ-9 Score"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

y_max <- max(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) - 8
p0.4 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = PTA, y = HL_PHQ9_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "PTA (dB HL)",
    y = "PHQ-9 Score",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HL_PHQ9_score, na.rm=TRUE) - 8
p0.5 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = HF.HL, y = HL_PHQ9_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +  
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "High frequency (dB HL)",
    y = "PHQ-9 Score"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = MoCA.Score, y = HL_PHQ9_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +   
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "MoCA Score",
       y = "PHQ-9 Score",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = thi_score, y = HL_PHQ9_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "THI Score",
       y = "PHQ-9 Score") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p2) <- names(dt)

p3 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = guef_score, y = HL_PHQ9_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "GUF Score",
       y = "PHQ-9 Score") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p3) <- names(dt)

p4 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = PTA, y = HL_PHQ9_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "PTA (dB HL)",
       y = "PHQ-9 Score",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p4) <- names(dt)

p5 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = HF.HL, y = HL_PHQ9_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "High frequency (dB HL)",
       y = "PHQ-9 Score") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p5) <- names(dt)

# pdf("figures/PHQ9_MOCA_THI_GUF_scatter_plot.pdf", width=25, height=21)
png("figures/PHQ9_MOCA_THI_GUF_scatter_plot.png", width=25*300, height=21*300, res=300)
plot_grid(p0.1, p0.2, p0.3, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p3[["gba_genes"]], p3[["synap_genes"]], p3[["BD_genes"]], p3[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]],
                        p3[["tinnitus_genes"]], p3[["tinnitus_synap_genes"]], p3[["hyperacusis_genes"]], p3[["hyperacusis_synap_genes"]],
                        nrow=7, ncol=4, align="hv", labels = c("A", "", "", "",
                                                          "B", "C", "D", "E",
                                                          "", "", "", "",
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", "", "",
                                                          "", "", "", "", ""))
dev.off()

# pdf("figures/PHQ9_HL_scatter_plot.pdf", width=25, height=15)
png("figures/PHQ9_HL_scatter_plot.png", width=25*300, height=15*300, res=300)
plot_grid(p0.4, p0.5, NULL, NULL,
                        p4[["gba_genes"]], p4[["synap_genes"]], p4[["BD_genes"]], p4[["PRKRA"]],
                        p5[["gba_genes"]], p5[["synap_genes"]], p5[["BD_genes"]], p5[["PRKRA"]],
                        p4[["tinnitus_genes"]], p4[["tinnitus_synap_genes"]], p4[["hyperacusis_genes"]], p4[["hyperacusis_synap_genes"]],
                        p5[["tinnitus_genes"]], p5[["tinnitus_synap_genes"]], p5[["hyperacusis_genes"]], p5[["hyperacusis_synap_genes"]],
                        nrow=5, ncol=4, align="hv", labels = c("A", "", "", "",
                                                          "B", "C", "D", "E",
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", "", ""))
dev.off()

# scatter plot for THI score by HL level
y_max <- max(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) - 30
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = thi_score, y = PTA, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "THI Score",
    y = "PTA (dB HL)",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) - 30
p0.2 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = thi_score, y = HF.HL, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "THI Score",
    y = "High requency (dB HL)"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = thi_score, y = PTA)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
   stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "THI Score",
       y = "PTA (dB HL)",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = thi_score, y = HF.HL)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
   stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "THI Score",
       y = "High frequency (dB HL)") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p2) <- names(dt)

# pdf("figures/HL_THI_scatter_plot_v2.pdf", width=25, height=15)
png("figures/HL_THI_scatter_plot_v2.png", width=25*300, height=15*300, res=300)
plot_grid(p0.1, p0.2, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]],
                        nrow=5, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", ""))
dev.off()

# scatter plot for GUF score by HL level
y_max <- max(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) - 20
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = guef_score, y = PTA, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "GUF Score",
    y = "PTA (dB HL)",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

y_max <- max(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$HF.HL, na.rm=TRUE) - 30
p0.2 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = guef_score, y = HF.HL, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "GUF Score",
    y = "High requency (dB HL)"
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = guef_score, y = PTA)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") + 
   stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "GUF Score",
       y = "PTA (dB HL)",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = guef_score, y = HF.HL)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "GUF Score",
       y = "High frequency (dB HL)") +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  )
})
names(p2) <- names(dt)

# pdf("figures/HL_GUF_scatter_plot_v2.pdf", width=25, height=15)
png("figures/HL_GUF_scatter_plot_v2.png", width=25*300, height=15*300, res=300)
plot_grid(p0.1, p0.2, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]],
                        nrow=5, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "", "", "", "",
                                                          "F", "G", "H", "I",
                                                          "", "", "", ""))
dev.off()

# scatter plot for PTA score by HFHL level
y_max <- max(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$PTA, na.rm=TRUE) - 20
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = HF.HL, y = PTA, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +  
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "High frequency (dB HL)",
    y = "PTA (dB HL)",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = HF.HL, y = PTA)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +   
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "High frequency (dB HL)",
       y = "PTA (dB HL)",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)


# pdf("figures/PTA_HFHL_scatter_plot.pdf", width=25, height=9)
png("figures/PTA_HFHL_scatter_plot.png", width=25*300, height=9*300, res=300)
plot_grid(p0.1, NULL, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        nrow=3, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "F", "G", "H", "I"))
dev.off()

# scatter plot for THI score by GUF level
y_max <- max(APOE_THI_GUEF_MOCA$guef_score, na.rm=TRUE) + 2
y_min <- min(APOE_THI_GUEF_MOCA$guef_score, na.rm=TRUE) - 10
p0.1 <- ggplot(APOE_THI_GUEF_MOCA[!is.na(APOE_THI_GUEF_MOCA$MoCA.level),], aes(x = thi_score, y = guef_score, color = MoCA.level)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  stat_poly_eq(
    aes(
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~"),
      color = MoCA.level,      # show text in same color as line
      group = MoCA.level
    ),
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4,
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(
    x = "THI Score",
    y = "GUF Score",
    title = "All samples (n = 294)",
  ) +
  theme_classic() +
  scale_color_manual(
    name = "MoCA status",   # Legend title
    values = c("Q1" = "#F8766D", "normal" = "#00BFC4"),  # colors for Q1 and normal
    labels = c("Q1" = "MCI (n = 75)", "normal" = "Non-MCI (n = 201)")  # legend labels
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  ggplot(x, aes(x = thi_score, y = guef_score)) +
  geom_point(alpha = 0.7, size=2, width=0.2, color = "#F8766D") +                        # scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#F8766D") +  
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",            # Or 0
    label.y = "bottom",
    size = 4, color = "#F8766D"
  ) +     
  labs(x = "THI Score",
       y = "GUF Score",
       title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme_classic() +
  scale_color_manual(values="#F8766D") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
})
names(p1) <- names(dt)


# pdf("figures/THI_GUF_scatter_plot.pdf", width=25, height=9)
png("figures/THI_GUF_scatter_plot.png", width=25*300, height=9*300, res=300)
plot_grid(p0.1, NULL, NULL, NULL,
                        p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]],
                        nrow=3, ncol=4, align="hv", labels = c("A", "", "", "", 
                                                          "B", "C", "D", "E", 
                                                          "F", "G", "H", "I"))
dev.off()

# boxplot for MoCA score by APOE risk and tinnitus level
p0.1 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = APOE_risk_1, y = MoCA.Score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  scale_y_continuous(breaks = seq(20, 30, by = 2)) +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "MoCA Score", title = "All samples\n(n = 294)") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
print(wilcox.test(MoCA.Score ~ APOE_risk_1, data = APOE_THI_GUEF_MOCA))

p0.2 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = APOE_risk_1, y = thi_score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "THI Score") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

print(wilcox.test(thi_score ~ APOE_risk_1, data = APOE_THI_GUEF_MOCA))

p0.3 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = APOE_risk_1, y = guef_score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "GUF Score") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
print(wilcox.test(guef_score ~ APOE_risk_1, data = APOE_THI_GUEF_MOCA))

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  print(n)
  print(wilcox.test(MoCA.Score ~ APOE_risk_1, data = x))
  ggplot(x, aes(x = APOE_risk_1, y = MoCA.Score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width=0.2, height=0, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "MoCA Score", title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
})
names(p1) <- names(dt)

p1[["PRKRA"]] <- p1[["PRKRA"]] + 
  scale_y_continuous(breaks = seq(23, 25, by = 1))

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  print(n)
  print(wilcox.test(thi_score ~ APOE_risk_1, data = x))
  ggplot(x, aes(x = APOE_risk_1, y = thi_score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "THI Score") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
})
names(p2) <- names(dt)

p3 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  print(n)
  print(wilcox.test(guef_score ~ APOE_risk_1, data = x))
  ggplot(x, aes(x = APOE_risk_1, y = guef_score, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "GUF Score") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
})
names(p3) <- names(dt)

pdf("figures/THI_GUEF_MOCA_box_plot_v3.pdf", width=13, height=15)
# png("figures/THI_GUEF_MOCA_box_plot_v3.png", width=13*300, height=15*300, res=300)
plot_grid(p0.1, p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p0.2, p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p0.3, p3[["gba_genes"]], p3[["synap_genes"]], p3[["BD_genes"]], p3[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]], NULL,
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]], NULL,
                        p3[["tinnitus_genes"]], p3[["tinnitus_synap_genes"]], p3[["hyperacusis_genes"]], p3[["hyperacusis_synap_genes"]], NULL,
                        nrow=6, ncol=5, align="hv", labels = c("A", "B", "C", "D", "E",
                                                          "", "", "", "", "",
                                                          "", "", "", "", "",
                                                          "F", "G", "H", "I", "",
                                                          "", "", "", "", "",
                                                          "", "", "", "", ""))
dev.off()

# # pdf("figures/THI_GUEF_MOCA_box_plot.pdf", width=13, height=15)
# png("figures/THI_GUEF_MOCA_box_plot.png", width=13*300, height=15*300, res=300)                                                      
# plot_grid(p0.1, p0.2, p0.3,
#                         p1[["gba_genes"]], p2[["gba_genes"]], p3[["gba_genes"]], 
#                         p1[["AD_genes"]], p2[["AD_genes"]], p3[["AD_genes"]],
#                         p1[["tinnitus_genes"]], p2[["tinnitus_genes"]], p3[["tinnitus_genes"]],
#                         p1[["tinnitus_synap_genes"]], p2[["tinnitus_synap_genes"]], p3[["tinnitus_synap_genes"]],
#                         nrow=5, ncol=3, align="hv", labels = c("A", "", "", 
#                                                           "B", "", "",
#                                                           "C", "", "",
#                                                           "D", "", "",
#                                                           "E", "", "",
#                                                           "F", "", "", 
#                                                           "G", "", ""))
# dev.off()

# boxplot for MoCA score by APOE risk and HL level
p0.1 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = APOE_risk_1, y = PTA, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "PTA (dB HL)", title = "All samples (n = 294)") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
print(wilcox.test(PTA ~ APOE_risk_1, data = APOE_THI_GUEF_MOCA))

p0.2 <- ggplot(APOE_THI_GUEF_MOCA, aes(x = APOE_risk_1, y = HF.HL, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "High frequency (dB HL)") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

print(wilcox.test(HF.HL ~ APOE_risk_1, data = APOE_THI_GUEF_MOCA))

p1 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  print(n)
  print(wilcox.test(PTA ~ APOE_risk_1, data = x))
  ggplot(x, aes(x = APOE_risk_1, y = PTA, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width=0.2, height=0, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "PTA (dB HL)", title = paste0(gene_types[[n]], "\n(# samples = ", nrow(x), ")")) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
})
names(p1) <- names(dt)

p2 <- lapply(names(dt), function(n) {
  x <- dt[[n]]
  print(n)
  print(wilcox.test(HF.HL ~ APOE_risk_1, data = x))
  ggplot(x, aes(x = APOE_risk_1, y = HF.HL, fill=APOE_risk_1)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(
    values = c("APOE E4" = "#F8766D", "Non-APOE E4" = "#00BFC4")
  ) +
  labs(x = "Group", y = "High frequency (dB HL)") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
})
names(p2) <- names(dt)

pdf("figures/HL_APOE_box_plot_v2.pdf", width=13, height=10)
# png("figures/HL_APOE_box_plot_v2.png", width=13*300, height=10*300, res=300)
plot_grid(p0.1, p1[["gba_genes"]], p1[["synap_genes"]], p1[["BD_genes"]], p1[["PRKRA"]],
                        p0.2, p2[["gba_genes"]], p2[["synap_genes"]], p2[["BD_genes"]], p2[["PRKRA"]],
                        p1[["tinnitus_genes"]], p1[["tinnitus_synap_genes"]], p1[["hyperacusis_genes"]], p1[["hyperacusis_synap_genes"]], NULL,
                        p2[["tinnitus_genes"]], p2[["tinnitus_synap_genes"]], p2[["hyperacusis_genes"]], p2[["hyperacusis_synap_genes"]], NULL,
                        nrow=4, ncol=5, align="hv", labels = c("A", "B", "C", "D", "E",
                                                          "", "", "", "", "",
                                                          "F", "G", "H", "I", "",
                                                          "", "", "", "", ""))
dev.off()

