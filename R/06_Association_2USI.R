#############################################################

# ASSOCIATION ANALYSIS WITHIN LINKAGE REGIONS 

#############################################################

# Working directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
       setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
     }

# Libraries
library(data.table)
library(dplyr)
library(qqman)

#############################################################

###                  GENOTYPE DATA 

#############################################################

# Create new bed, bim and fam files for association study
# Remove unrelated individuals

call1 <- paste0("plink1 --noweb --bfile ./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort",
              " --remove ./data/indiv_remove.txt --make-bed --out association_WISARD/KIN_assoc")
try(system(call1))

#############################################################

###                 COVARIATES

#############################################################


# Create the covariate file for GEMMA
pheno_data <- read.table("./data/2019-020_data.txt", header = TRUE, sep = "")
fam_data_original <- read.table("association_WISARD/KIN_assoc.fam", sep = "") 

cov <- dplyr::select(pheno_data, new_id, 
                            t722_smoking_habits_3m_BL,
                            t37_sample_age_BL,
                            t276_BMI_BL)
setnames(cov, old = c("new_id", 
                      "t722_smoking_habits_3m_BL",
                      "t37_sample_age_BL",
                      "t276_BMI_BL"),
         new = c("V2","SmokStat","Age","BMI"))
head(cov)

# Pre-processing

#### SmokStat
# Set missing as never
cov$SmokStat <- gsub("9999999", "1", cov$SmokStat)
# Set everything to never, except current smokers
cov$SmokStat <- gsub("2", "1", cov$SmokStat)
cov$SmokStat <- gsub("3", "1", cov$SmokStat)
# Now set 4 to 2
cov$SmokStat <- gsub("4", "2", cov$SmokStat)

# We need to match cov with fam_data
fam_cov <- merge(fam_data_original, cov, by = "V2", all.x = TRUE, 
                 sort = FALSE)[, union(names(fam_data_original), names(cov))]
head(fam_cov)

# Remove father and mother id
fam_cov$V3 <- NULL
fam_cov$V4 <- NULL

# Set names
setnames(fam_cov,old = c("V1","V2", "V5","V6"),
         new = c("FID", "IID", "Sex", "Aff_stat"))
# Check NA
table(is.na(fam_cov$BMI))

# Set missing as -9
fam_cov$SmokStat[is.na(fam_cov$SmokStat)] <- -9

### Age
# Set missing as -9
fam_cov$Age[is.na(fam_cov$Age)] <- -9

### BMI
# Set missing as -9
fam_cov$BMI[is.na(fam_cov$BMI)] <- -9
fam_cov$BMI[fam_cov$BMI == 7777777] <- -9

#############################################################

###                 QUANTITATIVE TRAITS

#############################################################

# Get quantitative phenotype from previous ped file
quant_trait <- fread("NEW_KIN_CHR_22.ped", data.table = F)
dim(quant_trait) # 1211 14834

# Select the microbiome part
quant_data <- quant_trait[,c(2,14809:14830)]
# Set taxa names
names(quant_data) <- c("V2","Clostridium_XlVa", "Roseburia", 
                       "Anaerostipes", "Escherichia/Shigella", 
                       "Dialister", "Alistipes", "Bacteroides", 
                       "Barnesiella", "Odoribacter", "Parabacteroides", 
                       "Bilophila", "Parasutterella", "Anaerotruncus", 
                       "Faecalibacterium", "Subdoligranulum", 
                       "Ruminococcus", "Clostridium_IV", "Flavonifractor", 
                       "Pseudoflavonifractor", "Oscillibacter", 
                       "Blautia", "Ruminococcus2")
# Set missing to -9
quant_data[quant_data == "X"] <- -9

# Merge microbiome with covariates
fam_cov$V2 <- gsub("2019-020_pre_", "",fam_cov$IID)

# Check if the two dfs are identical
identical(fam_cov$V2, quant_data$V2) # FALSE
identical(sort(fam_cov$V2), sort(quant_data$V2)) #TRUE
# So, the orderings are different. So, we merge and not cbind

# Now we can merge on the V2 variable
fam_cov_micr <- merge(fam_cov, quant_data, by= "V2", all.x = T, 
                      sort = F)[, union(names(fam_cov), names(quant_data))]
table(is.na(fam_cov_micr)) # FALSE all

# Write to file
write.table(fam_cov_micr, "association_WISARD/KIN_assoc_phen.txt", 
            sep = " ", quote = F, col.names = T, row.names = F)

# Perform GEMMA with covariate adjustment (Roseburia)

call3 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Roseburia --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/roseburia_gemma_cov")

try(system(call3))

# Results

GemmaRes_roseburia <- fread("./association_WISARD/roseburia_gemma_cov.gemma.res", 
                     sep = "\t", data.table = F, header = T)

# Manhattan Plot of Roseburia Association Test
assoc_res1 <- subset(GemmaRes_roseburia, !CHR == "X")
assoc_res2 <- subset(assoc_res1, !CHR == "Y")
assoc_res3 <- subset(assoc_res2, !CHR == "Un")
assoc_res4 <- subset(assoc_res3, !CHR == "XY")
assoc_res <- dplyr::select(assoc_res4, CHR, POS, VARIANT, P_SCORE)

assoc_res$CHR <- as.numeric(assoc_res$CHR)
table(is.na(assoc_res))
assoc_res_final <- na.omit(assoc_res)
manhattan(assoc_res_final, chr = "CHR", bp = "POS", snp = "VARIANT", p = "P_SCORE")

#####################

# Perform GEMMA with covariate adjustment (Closttidium_XIVa)

call4 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Clostridium_XlVa --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Clostridium_XlVa_gemma_cov")

try(system(call4))

#####################

# Perform GEMMA with covariate adjustment (Ruminococcus)

call5 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Ruminococcus --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Ruminococcus_gemma_cov")

try(system(call5))

#####################

# Perform GEMMA with covariate adjustment (Barnesiella)

call6 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Barnesiella --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Barnesiella_gemma_cov")

try(system(call6))

#####################

# Perform GEMMA with covariate adjustment (Pseudoflavonifractor)

call7 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Pseudoflavonifractor --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Pseudoflavonifractor_gemma_cov")

try(system(call7))

#####################

# Perform GEMMA with covariate adjustment (Parasutterella)

call8 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Parasutterella --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Parasutterella_gemma_cov")

try(system(call8))

#####################

# Perform GEMMA with covariate adjustment (Odoribacter)

call9 <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen.txt",
                " --pname Odoribacter --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out ./association_WISARD/Odoribacter_gemma_cov")

try(system(call9))


###############################################################################
##### ASSOCIATION ANALYSIS FOR SHANNON (ALPHA DIVERSITY)
###############################################################################

### QUANTITATIVE TRAITS TO GET SHANON 
# I#ve kept the variable name same as of now. Change it later
# Write a function
# Get quantitative phenotype from previous ped file
quant_trait <- fread("AD_KIN_CHR_22.ped", data.table = F)
dim(quant_trait) # 1211 14814

# Select the alpha diversity part
quant_data <- quant_trait[,c(2,14809:14810)]
# Set taxa names
names(quant_data) <- c("V2","Chao1","Shannon")
# Set missing to -9
quant_data[quant_data == "X"] <- -9

# Merge microbiome with covariates
fam_cov$V2 <- gsub("2019-020_pre_", "",fam_cov$IID)

# Check if the two dfs are identical
identical(fam_cov$V2, quant_data$V2) # FALSE
identical(sort(fam_cov$V2), sort(quant_data$V2)) #TRUE
# So, the orderings are different. So, we merge and not cbind

# Now we can merge on the V2 variable
fam_cov_micr <- merge(fam_cov, quant_data, by= "V2", all.x = T, 
                      sort = F)[, union(names(fam_cov), names(quant_data))]
table(is.na(fam_cov_micr)) # FALSE all

# Write to file
write.table(fam_cov_micr, "association_WISARD/KIN_assoc_phen_AD.txt", 
            sep = " ", quote = F, col.names = T, row.names = F)

# Perform GEMMA with covariate adjustment (Roseburia)

call_AD <- paste0("./wisard --bed association_WISARD/KIN_assoc.bed --sampvar association_WISARD/KIN_assoc_phen_AD.txt",
                " --pname Shannon --cname Sex,Aff_stat,SmokStat,Age,BMI --gemma --out association_WISARD/Shannon_gemma_cov")

try(system(call_AD))


###############################################################################

### EXTRACT LINKAGE REGION OF INTEREST (LROI)

###############################################################################

### Roseburia_CHR19

GemmaRes_roseburia <- data.table::fread("./association_WISARD/roseburia_gemma_cov.gemma.res", 
                            sep = "\t", data.table = F, header = T)
GemmaRes_roseburia <- subset(GemmaRes_roseburia, CHR == 19)

LROI_roseburia <- read.table("./Two-unit_support_interval/roseburia_2USI.txt", sep = "\t", header = T)

GEMMA_LROI <- GemmaRes_roseburia[GemmaRes_roseburia$POS %in% LROI_roseburia$BP,]
GEMMA_LROI <- dplyr::select(GEMMA_LROI,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI, "./Two-unit_support_interval/assoc_pvals/roseburia_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Clostridium_XlVa_CHR4 

GemmaRes_clost_4 <- data.table::fread("./association_WISARD/Clostridium_XlVa_gemma_cov.gemma.res", 
                                        sep = "\t", data.table = F, header = T)
GemmaRes_clost_4 <- subset(GemmaRes_clost_4, CHR == 4)

LROI_clost_4 <- read.table("./Two-unit_support_interval/clostridiumCHR4_2USI.txt", sep = "\t", header = T)

GEMMA_LROI1 <- GemmaRes_clost_4[GemmaRes_clost_4$POS %in% LROI_clost_4$BP,]
GEMMA_LROI1 <- dplyr::select(GEMMA_LROI1,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI1, "./Two-unit_support_interval/assoc_pvals/clostridium_CHR4_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Clostridium_XlVa_CHR14 (Two regions)

GemmaRes_clost_14 <- data.table::fread("./association_WISARD/Clostridium_XlVa_gemma_cov.gemma.res", 
                                      sep = "\t", data.table = F, header = T)
GemmaRes_clost_14 <- subset(GemmaRes_clost_14, CHR == 14)

# REGION 1
LROI_clost_14_1 <- read.table("./Two-unit_support_interval/clostridiumCHR14_2USI_region1.txt", sep = "\t", header = T)

GEMMA_LROI2_a <- GemmaRes_clost_14[GemmaRes_clost_14$POS %in% LROI_clost_14_1$BP,]
GEMMA_LROI2_a <- dplyr::select(GEMMA_LROI2_a,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI2_a, "./Two-unit_support_interval/assoc_pvals/clostridium_CHR14_gemma_2USI_region1.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# REGION 2
LROI_clost_14_2 <- read.table("./Two-unit_support_interval/clostridiumCHR14_2USI_region2.txt", sep = "\t", header = T)

GEMMA_LROI2_b <- GemmaRes_clost_14[GemmaRes_clost_14$POS %in% LROI_clost_14_2$BP,]
GEMMA_LROI2_b <- dplyr::select(GEMMA_LROI2_b,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI2_b, "./Two-unit_support_interval/assoc_pvals/clostridium_CHR14_gemma_2USI_region2.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Barnesiella_CHR4

GemmaRes_barne_4 <- data.table::fread("./association_WISARD/Barnesiella_gemma_cov.gemma.res", 
                                      sep = "\t", data.table = F, header = T)
GemmaRes_barne_4 <- subset(GemmaRes_barne_4, CHR == 4)

LROI_barne_4 <- read.table("./Two-unit_support_interval/barnesiellaCHR4_2USI.txt", sep = "\t", header = T)

GEMMA_LROI3 <- GemmaRes_barne_4[GemmaRes_barne_4$POS %in% LROI_barne_4$BP,]
GEMMA_LROI3 <- dplyr::select(GEMMA_LROI3,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI3, "./Two-unit_support_interval/assoc_pvals/barnesiella_CHR4_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Parasutterella_CHR14

GemmaRes_parasut_14 <- data.table::fread("./association_WISARD/Parasutterella_gemma_cov.gemma.res", 
                                       sep = "\t", data.table = F, header = T)
GemmaRes_parasut_14 <- subset(GemmaRes_parasut_14, CHR == 14)

LROI_parasut_14 <- read.table("./Two-unit_support_interval/parasutterellaCHR14_2USI.txt", sep = "\t", header = T)

GEMMA_LROI4 <- GemmaRes_parasut_14[GemmaRes_parasut_14$POS %in% LROI_parasut_14$BP,]
GEMMA_LROI4 <- dplyr::select(GEMMA_LROI4,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI4, "./Two-unit_support_interval/assoc_pvals/parasutterella_CHR14_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Pseudoflavonifractor_CHR7

GemmaRes_pseudo_7 <- data.table::fread("./association_WISARD/Pseudoflavonifractor_gemma_cov.gemma.res", 
                                         sep = "\t", data.table = F, header = T)
GemmaRes_pseudo_7 <- subset(GemmaRes_pseudo_7, CHR == 7)

LROI_pseudo_7 <- read.table("./Two-unit_support_interval/pseudoflavonifractorCHR7_2USI.txt", sep = "\t", header = T)

GEMMA_LROI5 <- GemmaRes_pseudo_7[GemmaRes_pseudo_7$POS %in% LROI_pseudo_7$BP,]
GEMMA_LROI5 <- dplyr::select(GEMMA_LROI5,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI5, "./Two-unit_support_interval/assoc_pvals/pseudoflavonifractor_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Ruminococcus_CHR16

GemmaRes_rumino_16 <- data.table::fread("./association_WISARD/Ruminococcus_gemma_cov.gemma.res", 
                                       sep = "\t", data.table = F, header = T)
GemmaRes_rumino_16 <- subset(GemmaRes_rumino_16, CHR == 16)

LROI_rumino_16 <- read.table("./Two-unit_support_interval/ruminococcusCHR16_2USI_region1&2.txt", sep = "\t", header = T)

GEMMA_LROI6 <- GemmaRes_rumino_16[GemmaRes_rumino_16$POS %in% LROI_rumino_16$BP,]
GEMMA_LROI6 <- dplyr::select(GEMMA_LROI6,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI6, "./Two-unit_support_interval/assoc_pvals/ruminococcus_gemma_2USI_region1&2.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Odoribacter_CHR22 (Region 1)

GemmaRes_odori1_22 <- data.table::fread("./association_WISARD/Odoribacter_gemma_cov.gemma.res", 
                                        sep = "\t", data.table = F, header = T)
GemmaRes_odori1_22 <- subset(GemmaRes_odori1_22, CHR == 22)

LROI_odori1_22 <- read.table("./Two-unit_support_interval/odoribacterCHR22_2USI_region1.txt", sep = "\t", header = T)

GEMMA_LROI7 <- GemmaRes_odori1_22[GemmaRes_odori1_22$POS %in% LROI_odori1_22$BP,]
GEMMA_LROI7 <- dplyr::select(GEMMA_LROI7,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI7, "./Two-unit_support_interval/assoc_pvals/odoribacter_gemma_2USI_region1.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Odoribacter_CHR22 (Region 2)

GemmaRes_odori2_22 <- data.table::fread("./association_WISARD/Odoribacter_gemma_cov.gemma.res", 
                                        sep = "\t", data.table = F, header = T)
GemmaRes_odori2_22 <- subset(GemmaRes_odori2_22, CHR == 22)

LROI_odori2_22 <- read.table("./Two-unit_support_interval/odoribacterCHR22_2USI_region2.txt", sep = "\t", header = T)

GEMMA_LROI8 <- GemmaRes_odori2_22[GemmaRes_odori2_22$POS %in% LROI_odori2_22$BP,]
GEMMA_LROI8 <- dplyr::select(GEMMA_LROI8,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI8, "./Two-unit_support_interval/assoc_pvals/odoribacter_gemma1_2USI_region2.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# REGION 3

LROI_odori3_22 <- read.table("./Two-unit_support_interval/odoribacterCHR22_2USI_region3.txt", sep = "\t", header = T)
GEMMA_LROI9 <- GemmaRes_odori2_22[GemmaRes_odori2_22$POS %in% LROI_odori3_22$BP,]
GEMMA_LROI9 <- dplyr::select(GEMMA_LROI9,CHR,VARIANT,POS,ALT,P_LRT)

# write to file
write.table(GEMMA_LROI9, "./Two-unit_support_interval/assoc_pvals/odoribacter_gemma1_2USI_region3.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Shannon_Index_CHR3

GemmaRes_shannon_3 <- data.table::fread("association_WISARD/Shannon_gemma_cov.gemma.res", 
                                         sep = "\t", data.table = F, header = T)
GemmaRes_shannon_3 <- subset(GemmaRes_shannon_3, CHR == 3)

LROI_shannon_3 <- read.table("./Two-unit_support_interval/shannonCHR3_2USI.txt", sep = "\t", header = T)

GEMMA_LROI_10 <- GemmaRes_shannon_3[GemmaRes_shannon_3$POS %in% LROI_shannon_3$BP,]
GEMMA_LROI_10 <- dplyr::select(GEMMA_LROI_10,CHR,VARIANT,POS,ALT,P_LRT)
# write to file
write.table(GEMMA_LROI_10, "./Two-unit_support_interval/assoc_pvals/shannon_gemma_2USI.txt", sep = "\t", quote = F, row.names = F, col.names = T)





