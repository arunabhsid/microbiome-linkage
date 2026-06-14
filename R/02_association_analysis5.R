
#############################################################

### PAIRWISE GENETIC & MICROBIOME SIMILARITY ANALYSIS

#############################################################

####   ALL GENERA FROM KINDRED

# Run this script from the folder that contains it (the folder with the ./data subfolder).

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


# Libraries
library(phyloseq)
library(reshape2)
library(ggplot2)
library(dplyr)

all_genera_KIN <- readRDS("./data/KIN_genera_all.RDS")


# Have a look at the resulting dataset
ASV_to_genus <- as.data.frame(tax_table(all_genera_KIN)[,6])
dim(ASV_to_genus)

# Calculate distance
KIN_genera_dist <- distance(all_genera_KIN, method = "bray")
KIN_genera_dist <- as.matrix(KIN_genera_dist)
KIN_genera_dist[1:5,1:5]

# Convert it to a tabular format

df_KIN_genera_dist <- melt(as.matrix(KIN_genera_dist), varnames = c("row", "col"))
head(df_KIN_genera_dist)
# Concatenate row and col from df
df_KIN_genera_dist$PAIR <- paste0(df_KIN_genera_dist$row,"-",df_KIN_genera_dist$col)
head(df_KIN_genera_dist)
# Remove unwanted columns
bray_dist_KIN_genera <- df_KIN_genera_dist[,c(4,3)]
head(bray_dist_KIN_genera)


##########################

# GENETIC DATA

##########################

std_ibd_WIDE <- readRDS(file = "./data/std_IBD_WIDE.Rds")

##########################

# COVARIATES

##########################

load(file = "./data/Covariates_pairwise.RData")


##########################

# MERGE

##########################
### Merging genetic, microbiome data and covariates

# 1_Microbiome
data_selectGenera <- merge(std_ibd_WIDE, bray_dist_KIN_genera, by= "PAIR")
dim(data_selectGenera)   # 364 191
# 2_Sex
data_selectGenera2 <- merge(data_selectGenera, sex_df_final, by= "PAIR")
dim(data_selectGenera2)   # 364 192
# 3_Disease Status
data_selectGenera3 <- merge(data_selectGenera2, disStat_df_final, by= "PAIR")
dim(data_selectGenera3)   # 364 193
# 4_Smoking status
data_selectGenera4 <- merge(data_selectGenera3, SmoStat_df_final, by= "PAIR")
dim(data_selectGenera4)   # 364 194
# 5_Age difference
data_selectGenera5 <- merge(data_selectGenera4, Age_df_final, by= "PAIR")
dim(data_selectGenera5)   # 364 195
# 6_BMI
data_selectGenera6 <- merge(data_selectGenera5, BMI_df_final, by= "PAIR")
dim(data_selectGenera6)   # 364 196
# 7_Household (Geo coordinate)
data_BC_selectGenera <- merge(data_selectGenera6, GEO_dist, by= "PAIR")
dim(data_BC_selectGenera)   # 364 197

# Remove NA rows
data_BC_selectGenera <- data_BC_selectGenera[!grepl("NA", data_BC_selectGenera$DISSTAT),]
dim(data_BC_selectGenera)   # 3 rows removed
names(data_BC_selectGenera)[189:197]
dim(data_BC_selectGenera) # 361 197

# DISCORDANT PAIRS
data_BC_selectGenera_Disc <- dplyr::filter(data_BC_selectGenera, DISSTAT == "unaff-aff")
dim(data_BC_selectGenera_Disc) # 187 197

# CONCORDANT PAIRS (ALL)
data_BC_selectGenera_Conc <- dplyr::filter(data_BC_selectGenera, DISSTAT == "aff-aff" | DISSTAT == "unaff-unaff")
dim(data_BC_selectGenera_Conc) # 174 197

# CONCORDANT PAIRS (AFFECTED)
data_BC_selectGenera_Conc_aff <- dplyr::filter(data_BC_selectGenera, DISSTAT == "aff-aff")
dim(data_BC_selectGenera_Conc_aff) # 16 197

# CONCORDANT PAIRS (UNAFFECTED)
data_BC_selectGenera_Conc_unaff <- dplyr::filter(data_BC_selectGenera, DISSTAT == "unaff-unaff")
dim(data_BC_selectGenera_Conc_unaff) # 158 197

# Subset the SNPs and covariates (ALL PAIRS)
snps_selectGenera <- data_BC_selectGenera[,c(2:190)]
dat_selectGenera <- data_BC_selectGenera[,c(1,191:197)]
data.table::setnames(dat_selectGenera, old = c("value.y","value.x"), new = c("household","betadiv"))

# Subset the SNPs and covariates (DISCORDANT PAIRS)
snps_selectGenera_disc <- data_BC_selectGenera_Disc[,c(2:190)]
dat_selectGenera_disc <- data_BC_selectGenera_Disc[,c(1,191:197)]
data.table::setnames(dat_selectGenera_disc, old = c("value.y","value.x"), new = c("household","betadiv"))

# Subset the SNPs and covariates (ALL CONCORDANT PAIRS)
snps_selectGenera_conc <- data_BC_selectGenera_Conc[,c(2:190)]
dat_selectGenera_conc <- data_BC_selectGenera_Conc[,c(1,191:197)]
data.table::setnames(dat_selectGenera_conc, old = c("value.y","value.x"), new = c("household","betadiv"))

# Subset the SNPs and covariates (CONCORDANT UNAFFECTED)
snps_selectGenera_conc_unaff <- data_BC_selectGenera_Conc_unaff[,c(2:190)]
dat_selectGenera_conc_unaff <- data_BC_selectGenera_Conc_unaff[,c(1,191:197)]
data.table::setnames(dat_selectGenera_conc_unaff, old = c("value.y","value.x"), new = c("household","betadiv"))

# Subset the SNPs and covariates (CONCORDANT AFFECTED)
snps_selectGenera_conc_aff <- data_BC_selectGenera_Conc_aff[,c(2:190)]
dat_selectGenera_conc_aff <- data_BC_selectGenera_Conc_aff[,c(1,191:197)]
data.table::setnames(dat_selectGenera_conc_aff, old = c("value.y","value.x"), new = c("household","betadiv"))

##############################################################

###         LINEAR MODEL (ALL PAIRS) - SKIP

##############################################################

pvalue_snp_selectGenera <- numeric(ncol(snps_selectGenera))
names(pvalue_snp_selectGenera) <- colnames(snps_selectGenera)

for (i in 1:ncol(snps_selectGenera)) {
  fit_selectGenera <- lm(dat_selectGenera$betadiv ~ snps_selectGenera[,i]
                         + factor(dat_selectGenera$DISSTAT) 
                         + factor(dat_selectGenera$SMOSTAT)
                         + dat_selectGenera$AGE_DIFF
                         + dat_selectGenera$BMI_DIFF
                         + factor(dat_selectGenera$SEX)
                         + factor(dat_selectGenera$household))
  pvalue_snp_selectGenera[i] <- anova(fit_selectGenera)['snps_selectGenera[, i]','Pr(>F)']
  
  
}

sorted_pvalues_selectGenera <- data.frame(pvalue = sort(pvalue_snp_selectGenera)[1:20])
print(sorted_pvalues_selectGenera) 

# Adjusting for multiple testing
qvalue_selectGenera <- p.adjust(pvalue_snp_selectGenera, 'fdr')
sorted_qvalues_selectGenera <- data.frame(qvalue = sort(qvalue_selectGenera)[1:20])
print(sorted_qvalues_selectGenera)

##############################################################

###         LINEAR MODEL (DISCORDANT PAIRS)

##############################################################

pvalue_snp_selectGenera_disc <- numeric(ncol(snps_selectGenera_disc))
names(pvalue_snp_selectGenera_disc) <- colnames(snps_selectGenera_disc)

for (i in 1:ncol(snps_selectGenera_disc)) {
  fit_selectGenera_disc <- lm(dat_selectGenera_disc$betadiv ~ snps_selectGenera_disc[,i]
                              + factor(dat_selectGenera_disc$SMOSTAT)
                              + dat_selectGenera_disc$AGE_DIFF
                              + factor(dat_selectGenera_disc$SEX)
                              + dat_selectGenera_disc$BMI_DIFF
                              + dat_selectGenera_disc$household)
  pvalue_snp_selectGenera_disc[i] <- anova(fit_selectGenera_disc)['snps_selectGenera_disc[, i]','Pr(>F)']
}

sorted_pvalues_selectGenera_disc <- data.frame(pvalue = sort(pvalue_snp_selectGenera_disc)[1:20])
print(sorted_pvalues_selectGenera_disc) 

# Adjusting for multiple testing
qvalue_selectGenera_disc <- p.adjust(pvalue_snp_selectGenera_disc, 'fdr')
sorted_qvalues_selectGenera_disc <- data.frame(qvalue = sort(qvalue_selectGenera_disc)[1:20])
print(sorted_qvalues_selectGenera_disc)

# Table for paper
all_snps <- read.table("./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort.bim")
all_snps <- dplyr::select(all_snps, V2,V1,V4)
all_snps$CHR_POS <- paste(all_snps$V1,all_snps$V4, sep= ":")
data.table::setnames(all_snps, old = "V2", new = "SNP")
all_snps <- dplyr::select(all_snps, SNP, CHR_POS)
all_snps$SNP <- gsub("GSA-","", all_snps$SNP)
supp_tab_3 <- data.frame(SNP = rownames(sorted_pvalues_selectGenera_disc), 
                         p = sorted_pvalues_selectGenera_disc, 
                         q = sorted_qvalues_selectGenera_disc)
supp_tab_3 <- supp_tab_3 %>% mutate(across(where(is.numeric), ~ round(.,6)))
# merge
dim(supp_tab_3)
dim(all_snps)
suppl_table_3 <- merge(supp_tab_3, all_snps, by = "SNP", sort = F)
write.csv(suppl_table_3, file = "./results/Suppl_table_3.csv", quote = F, row.names = F)

# rs11741861: qvalue = 0.005053516
lm_rs11741861 <- lm(dat_selectGenera_disc$betadiv ~ snps_selectGenera_disc$rs11741861
                    + factor(dat_selectGenera_disc$SEX)
                    + factor(dat_selectGenera_disc$SMOSTAT)
                    + dat_selectGenera_disc$AGE_DIFF
                    + dat_selectGenera_disc$BMI_DIFF
                    + dat_selectGenera_disc$household)
summary(lm_rs11741861)
# Save output
sink("./results/lm_rs11741861.txt")
print(summary(lm_rs11741861))
print("----Confidence Interval----")
# Confidence Interval
confint(lm_rs11741861)
sink()

# Plot
head(data_BC_selectGenera_Disc)

data.table::setnames(data_BC_selectGenera_Disc,
                     old = "value.x", new= "betadiv")
# Subset 
Plot_df <- dplyr::select(data_BC_selectGenera_Disc,PAIR, 
                         rs11741861, betadiv, DISSTAT, SMOSTAT)

jpeg("./plots/rs11741861.jpg",height=6,width=4,units='in',res=600)
plot1 <- ggplot(Plot_df, aes(x=rs11741861, y=betadiv)) + 
  geom_point() + xlab("Excess IBD Sharing") +
  ylab("Beta Diversity (Bray Curtis)")+
  ggtitle("rs11741861") +
  geom_smooth(method=lm, color='#2C3E50') + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
plot1
dev.off()
##############################################################

###         LINEAR MODEL (CONCORDANT PAIRS)

##############################################################

pvalue_snp_selectGenera_conc <- numeric(ncol(snps_selectGenera_conc))
names(pvalue_snp_selectGenera_conc) <- colnames(snps_selectGenera_conc)

for (i in 1:ncol(snps_selectGenera_conc)) {
  fit_selectGenera_conc <- lm(dat_selectGenera_conc$betadiv ~ snps_selectGenera_conc[,i]
                              + factor(dat_selectGenera_conc$SMOSTAT)
                              + dat_selectGenera_conc$AGE_DIFF
                              + factor(dat_selectGenera_conc$SEX))
  pvalue_snp_selectGenera_conc[i] <- anova(fit_selectGenera_conc)['snps_selectGenera_conc[, i]','Pr(>F)']
  
  
}

sorted_pvalues_selectGenera_conc <- data.frame(pvalue = sort(pvalue_snp_selectGenera_conc)[1:20])
print(sorted_pvalues_selectGenera_conc) 

# Adjusting for multiple testing
qvalue_selectGenera_conc <- p.adjust(pvalue_snp_selectGenera_conc, 'fdr')
sorted_qvalues_selectGenera_conc <- data.frame(qvalue = sort(qvalue_selectGenera_conc)[1:20])
print(sorted_qvalues_selectGenera_conc)

# NO SIGNIFICANT SNPs

##############################################################

###         LINEAR MODEL (CONCORDANT UNAFFECTED)

##############################################################

pvalue_snp_selectGenera_conc_unaff <- numeric(ncol(snps_selectGenera_conc_unaff))
names(pvalue_snp_selectGenera_conc_unaff) <- colnames(snps_selectGenera_conc_unaff)

for (i in 1:ncol(snps_selectGenera_conc_unaff)) {
  fit_selectGenera_conc_unaff <- lm(dat_selectGenera_conc_unaff$betadiv ~ snps_selectGenera_conc_unaff[,i]
                              + factor(dat_selectGenera_conc_unaff$SMOSTAT)
                              + dat_selectGenera_conc_unaff$AGE_DIFF
                              + factor(dat_selectGenera_conc_unaff$SEX))
  pvalue_snp_selectGenera_conc_unaff[i] <- anova(fit_selectGenera_conc_unaff)['snps_selectGenera_conc_unaff[, i]','Pr(>F)']
  
  
}

sorted_pvalues_selectGenera_conc_unaff <- data.frame(pvalue = sort(pvalue_snp_selectGenera_conc_unaff)[1:20])
print(sorted_pvalues_selectGenera_conc_unaff) 

# Adjusting for multiple testing
qvalue_selectGenera_conc_unaff <- p.adjust(pvalue_snp_selectGenera_conc_unaff, 'fdr')
sorted_qvalues_selectGenera_conc_unaff <- data.frame(qvalue = sort(qvalue_selectGenera_conc_unaff)[1:20])
print(sorted_qvalues_selectGenera_conc_unaff)

# NO SIGNIFICANT SNPs

##############################################################

###         LINEAR MODEL (CONCORDANT AFFECTED)

##############################################################

pvalue_snp_selectGenera_conc_aff <- numeric(ncol(snps_selectGenera_conc_aff))
names(pvalue_snp_selectGenera_conc_aff) <- colnames(snps_selectGenera_conc_aff)

for (i in 1:ncol(snps_selectGenera_conc_aff)) {
  fit_selectGenera_conc_aff <- lm(dat_selectGenera_conc_aff$betadiv ~ snps_selectGenera_conc_aff[,i]
                                    + factor(dat_selectGenera_conc_aff$SMOSTAT)
                                    + dat_selectGenera_conc_aff$AGE_DIFF
                                    + factor(dat_selectGenera_conc_aff$SEX))
  pvalue_snp_selectGenera_conc_aff[i] <- anova(fit_selectGenera_conc_aff)['snps_selectGenera_conc_aff[, i]','Pr(>F)']
  
  
}

sorted_pvalues_selectGenera_conc_aff <- data.frame(pvalue = sort(pvalue_snp_selectGenera_conc_aff)[1:20])
print(sorted_pvalues_selectGenera_conc_aff) 

# Adjusting for multiple testing
qvalue_selectGenera_conc_aff <- p.adjust(pvalue_snp_selectGenera_conc_aff, 'fdr')
sorted_qvalues_selectGenera_conc_aff <- data.frame(qvalue = sort(qvalue_selectGenera_conc_aff)[1:20])
print(sorted_qvalues_selectGenera_conc_aff)

# NO SIGNIFICANT SNPs

# Session Info
writeLines(capture.output(sessionInfo()), "./results/sessionInfo.txt")
