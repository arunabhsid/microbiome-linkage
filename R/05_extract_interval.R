#############################################################

# EXTRACTING LINKAGE REGIONS (Two-unit support interval)

#############################################################

# Libraries
library(ggplot2)
library(grid)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
}

# NOTE ON WORKFLOW:
# For each trait/region this script writes the full LOD profile to
# ./results/two_unit_interval/whole_region/. The two-unit support intervals
# were then determined manually from those LOD curves (SNPs within two LOD
# units of each peak) and saved as <trait>_2USI.txt in./results/two_unit_interval/.
# This script reads those interval files back in to verify and plot them.

## 2 LOD DROP from the highest LOD score on either side
## Arrange according to SNP-basepair position

### ROSEBURIA - CHR 19

LOD_CHR_19 <- read.table("./results/linkage_markernames/merlin-regress-chr19.tbl", 
                         sep = "\t", header = T)
LOD_Roseburia_CHR_19 <- subset(LOD_CHR_19, PHENOTYPE == "Trait: Roseburia")
dim(LOD_Roseburia_CHR_19) # 10875     8
LOD_Roseburia_CHR_19$PHENOTYPE <- gsub("Trait:", "", LOD_Roseburia_CHR_19$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_19_MAP <- read.table("KIN_CHR_19.map", sep = "")
CHR_19_MAP$V1 <- NULL
data.table::setnames(CHR_19_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Roseburia_CHR_19) # 10875     8
dim(CHR_19_MAP) # 11114     3
ROS_CHR_region <- merge(LOD_Roseburia_CHR_19, CHR_19_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Roseburia_CHR19.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.781

# Read file
Roseburia_2USI <- read.table("./results/two_unit_interval/roseburia_2USI.txt", 
                             sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Roseburia_2USI$BP)
dim(Roseburia_2USI) # 698   6

# Plot the region
ggplot(data = Roseburia_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### CLOSTRIDIUM_XIVa - CHR 14

LOD_CHR_14 <- read.table("./results/linkage_markernames/merlin-regress-chr14.tbl", 
                         sep = "\t", header = T)
LOD_Clostridium_XlVa_CHR_14 <- subset(LOD_CHR_14, PHENOTYPE == "Trait: Clostridium_XlVa")
dim(LOD_Clostridium_XlVa_CHR_14) # 16195     8
LOD_Clostridium_XlVa_CHR_14$PHENOTYPE <- gsub("Trait:", "", LOD_Clostridium_XlVa_CHR_14$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_14_MAP <- read.table("KIN_CHR_14.map", sep = "")
CHR_14_MAP$V1 <- NULL
data.table::setnames(CHR_14_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Clostridium_XlVa_CHR_14) # 16195     8
dim(CHR_14_MAP) # 16585     3
CLOS_CHR_region <- merge(LOD_Clostridium_XlVa_CHR_14, CHR_14_MAP, by = "POS", all.x = T, sort = F)
CLOS_CHR_region2 <- CLOS_CHR_region[order(CLOS_CHR_region$BP),]
CLOS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = CLOS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()

# Write to file (for manual analysis)
write.table(CLOS_CHR_region2, "./results/two_unit_interval/whole_region/Clostridium_XlVa_CHR14.txt", sep = "\t", quote = F, col.names = T, row.names = F)

### COMMENT
# There will be multiple regions!

# Find all regions where LOD is greater than 3
CLOS_LOD3 <- dplyr::filter(CLOS_CHR_region2, LOD >= 3)
# Plot
ggplot(data = CLOS_LOD3, aes(x=BP,y=LOD)) + geom_line()

# Find the maximum LOD in one region
max(CLOS_CHR_region2$LOD) # 4.387

# Read file 1
Clostridium14_2USI <- read.table("./results/two_unit_interval/clostridiumCHR14_2USI_region1.txt", 
                             sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Clostridium14_2USI$BP)
dim(Clostridium14_2USI) # 1221   6

# Plot the region
ggplot(data = Clostridium14_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!

# Another peak: LOD 3.597

# Read file 2
Clostridium14_2USI_2 <- read.table("./results/two_unit_interval/clostridiumCHR14_2USI_region2.txt", 
                                 sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Clostridium14_2USI_2$BP)
dim(Clostridium14_2USI_2) # 2139   6

# Plot the region
ggplot(data = Clostridium14_2USI_2, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### CLOSTRIDIUM_XIVa - CHR 4

LOD_CHR_4 <- read.table("./results/linkage_markernames/merlin-regress-chr04.tbl", 
                        sep = "\t", header = T)
LOD_Clostridium_XlVa_CHR_4 <- subset(LOD_CHR_4, PHENOTYPE == "Trait: Clostridium_XlVa")
dim(LOD_Clostridium_XlVa_CHR_4) # 30186     8
LOD_Clostridium_XlVa_CHR_4$PHENOTYPE <- gsub("Trait:", "", LOD_Clostridium_XlVa_CHR_4$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_4_MAP <- read.table("KIN_CHR_4.map", sep = "")
CHR_4_MAP$V1 <- NULL
data.table::setnames(CHR_4_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Clostridium_XlVa_CHR_4) # 30186     8
dim(CHR_4_MAP) # 32354     3
ROS_CHR_region <- merge(LOD_Clostridium_XlVa_CHR_4, CHR_4_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line() # One peak
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.022

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Clostridium_XlVa_CHR4.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Read file 
Clostridium4_2USI <- read.table("./results/two_unit_interval/clostridiumCHR4_2USI.txt", 
                                   sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Clostridium4_2USI$BP)
dim(Clostridium4_2USI) # 528   6

# Plot the region
ggplot(data = Clostridium4_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### Barnesiella - CHR 4

LOD_CHR_4 <- read.table("./results/linkage_markernames/merlin-regress-chr04.tbl", 
                        sep = "\t", header = T)
LOD_Barnesiella_CHR_4 <- subset(LOD_CHR_4, PHENOTYPE == "Trait: Barnesiella")
dim(LOD_Barnesiella_CHR_4) # 30186    8
LOD_Barnesiella_CHR_4$PHENOTYPE <- gsub("Trait:", "", LOD_Barnesiella_CHR_4$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_4_MAP <- read.table("KIN_CHR_4.map", sep = "")
CHR_4_MAP$V1 <- NULL
data.table::setnames(CHR_4_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Barnesiella_CHR_4) # 30186     8
dim(CHR_4_MAP) # 32354     3
ROS_CHR_region <- merge(LOD_Barnesiella_CHR_4, CHR_4_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.241

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Barnesiella_CHR4.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Read file 
Barnesiella4_2USI <- read.table("./results/two_unit_interval/barnesiellaCHR4_2USI.txt", 
                                sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Barnesiella4_2USI$BP)
dim(Barnesiella4_2USI) # 1821   6

# Plot the region
ggplot(data = Barnesiella4_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### Ruminococcus - CHR 16 (Two regions together)

LOD_CHR_16 <- read.table("./results/linkage_markernames/merlin-regress-chr16.tbl", 
                         sep = "\t", header = T)
LOD_Ruminococcus_CHR_16 <- subset(LOD_CHR_16, PHENOTYPE == "Trait: Ruminococcus")
dim(LOD_Ruminococcus_CHR_16) # 16442    8
LOD_Ruminococcus_CHR_16$PHENOTYPE <- gsub("Trait:", "", LOD_Ruminococcus_CHR_16$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_16_MAP <- read.table("KIN_CHR_16.map", sep = "")
CHR_16_MAP$V1 <- NULL
data.table::setnames(CHR_16_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Ruminococcus_CHR_16) # 16442     8
dim(CHR_16_MAP) # 16986     3
ROS_CHR_region <- merge(LOD_Ruminococcus_CHR_16, CHR_16_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 4.421

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Ruminococcus_CHR16.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Read file 
Ruminococcus_2USI <- read.table("./results/two_unit_interval/ruminococcusCHR16_2USI_region1&2.txt", 
                                sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Ruminococcus_2USI$BP)
dim(Ruminococcus_2USI) # 2391   6

# Plot the region
vertical_lines1 <- c(56674237,61100551)
ggplot(data = Ruminococcus_2USI, aes(x=BP,y=LOD)) + geom_line() + 
  geom_segment(aes(x=56674237, xend=61100551, y=1,yend=1),arrow=arrow(length=unit(0.3,"cm"),ends ="both"),size=0.5) + 
  geom_segment(aes(x=53851304, xend=65460917, y=0.9,yend=0.9),arrow=arrow(length=unit(0.3,"cm"),ends ="both"),size=0.5) +
  geom_segment(size=0.2,aes(x=58661456,xend=58661456,y=4.421,yend=1), linetype=2) +
  geom_segment(size=0.2,aes(x=61598734,xend=61598734,y=3.113,yend=0.9), linetype=2)

# Looks good!
###

### Pseudoflavonifractor - CHR 7 

LOD_CHR_7 <- read.table("./results/linkage_markernames/merlin-regress-chr07.tbl", 
                        sep = "\t", header = T)
LOD_Pseudoflavonifractor_CHR_7 <- subset(LOD_CHR_7, PHENOTYPE == "Trait: Pseudoflavonifractor")
dim(LOD_Pseudoflavonifractor_CHR_7) # 26719    8
LOD_Pseudoflavonifractor_CHR_7$PHENOTYPE <- gsub("Trait:", "", LOD_Pseudoflavonifractor_CHR_7$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_7_MAP <- read.table("KIN_CHR_7.map", sep = "")
CHR_7_MAP$V1 <- NULL
data.table::setnames(CHR_7_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Pseudoflavonifractor_CHR_7) # 26719     8
dim(CHR_7_MAP) # 28269     3
ROS_CHR_region <- merge(LOD_Pseudoflavonifractor_CHR_7, CHR_7_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.133

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Pseudoflavonifractor_CHR7.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Read file 
Pseudoflavonifractor_2USI <- read.table("./results/two_unit_interval/pseudoflavonifractorCHR7_2USI.txt", 
                                sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Pseudoflavonifractor_2USI$BP)
dim(Pseudoflavonifractor_2USI) # 1995   6

# Plot the region
ggplot(data = Pseudoflavonifractor_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### Parasutterella - CHR 14

LOD_CHR_14 <- read.table("./results/linkage_markernames/merlin-regress-chr14.tbl", 
                         sep = "\t", header = T)
LOD_Parasutterella_CHR_14 <- subset(LOD_CHR_14, PHENOTYPE == "Trait: Parasutterella")
dim(LOD_Parasutterella_CHR_14) # 16195    8
LOD_Parasutterella_CHR_14$PHENOTYPE <- gsub("Trait:", "", LOD_Parasutterella_CHR_14$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_14_MAP <- read.table("KIN_CHR_14.map", sep = "")
CHR_14_MAP$V1 <- NULL
data.table::setnames(CHR_14_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Parasutterella_CHR_14) # 16195     8
dim(CHR_14_MAP) # 16585     3
ROS_CHR_region <- merge(LOD_Parasutterella_CHR_14, CHR_14_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.018

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Parasutterella_CHR14.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Read file 
Parasutterella_2USI <- read.table("./results/two_unit_interval/parasutterellaCHR14_2USI.txt", 
                                        sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Parasutterella_2USI$BP)
dim(Parasutterella_2USI) # 986   6

# Plot the region
ggplot(data = Parasutterella_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### Odoribacter - CHR 22 (DISCUSS-- DO IT AFTER THE MEETINGs)

LOD_CHR_22 <- read.table("./results/linkage_markernames/merlin-regress-chr22.tbl", 
                         sep = "\t", header = T)
LOD_Odoribacter_CHR_22 <- subset(LOD_CHR_22, PHENOTYPE == "Trait: Odoribacter")
dim(LOD_Odoribacter_CHR_22) # 7313    8
LOD_Odoribacter_CHR_22$PHENOTYPE <- gsub("Trait:", "", LOD_Odoribacter_CHR_22$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_22_MAP <- read.table("KIN_CHR_22.map", sep = "")
CHR_22_MAP$V1 <- NULL
data.table::setnames(CHR_22_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Odoribacter_CHR_22) # 7313     8
dim(CHR_22_MAP) # 7401     3
ROS_CHR_region <- merge(LOD_Odoribacter_CHR_22, CHR_22_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()

# Write to file (for manual analysis)
write.table(ROS_CHR_region2, "./results/two_unit_interval/whole_region/Odoribacter_CHR22.txt", sep = "\t", quote = F, col.names = T, row.names = F)

### COMMENT
# There will be multiple regions!

# Find all regions where LOD is greater than 3
ODORI_LOD3 <- dplyr::filter(ROS_CHR_region2, LOD >= 3)
# Plot
ggplot(data = ODORI_LOD3, aes(x=BP,y=LOD)) + geom_line()
# The max of all LODs
max(ROS_CHR_region2$LOD) # 4.311

# Read file 1
Odoribacter_2USI_1 <- read.table("./results/two_unit_interval/odoribacterCHR22_2USI_region1.txt", 
                                  sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Odoribacter_2USI_1$BP)
dim(Odoribacter_2USI_1) # 22   6

# Plot the region
ggplot(data = Odoribacter_2USI_1, aes(x=BP,y=LOD)) + geom_line() # Looks good!

# Read file 2
Odoribacter_2USI_2 <- read.table("./results/two_unit_interval/odoribacterCHR22_2USI_region2.txt", 
                                 sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Odoribacter_2USI_2$BP)
dim(Odoribacter_2USI_2) # 763   6

# Plot the region
ggplot(data = Odoribacter_2USI_2, aes(x=BP,y=LOD)) + geom_line() # Looks good!

# Read file 3
Odoribacter_2USI_3 <- read.table("./results/two_unit_interval/odoribacterCHR22_2USI_region3.txt", 
                                 sep = "\t", header = T)
# Check if the region is not sorted
is.unsorted(Odoribacter_2USI_3$BP)
dim(Odoribacter_2USI_3) # 281   6

# Plot the region
ggplot(data = Odoribacter_2USI_3, aes(x=BP,y=LOD)) + geom_line() # Looks good!
###

### ALPHA DIVERSITY SHANNON INDEX - CHR 3

LOD_CHR_3 <- read.table("./results/linkage_AD/merlin-regress-chr03.tbl", 
                        sep = "\t", header = T)
LOD_Shannon_CHR_3 <- subset(LOD_CHR_3, PHENOTYPE == "Trait: Shannon")
dim(LOD_Shannon_CHR_3) # 32057     8
LOD_Shannon_CHR_3$PHENOTYPE <- gsub("Trait:", "", LOD_Shannon_CHR_3$PHENOTYPE)

# Map the SNPs to basepair position using the MAP fle
# Map file has both cM and base pair positions

CHR_3_MAP <- read.table("KIN_CHR_3.map", sep = "")
CHR_3_MAP$V1 <- NULL
data.table::setnames(CHR_3_MAP, old = c("V2","V3","V4"),
                     new = c("POS","cM","BP"))
# Merge 
dim(LOD_Shannon_CHR_3) # 32057     8
dim(CHR_3_MAP) # 34873     3
SHA_CHR_region <- merge(LOD_Shannon_CHR_3, CHR_3_MAP, by = "POS", all.x = T, sort = F)
SHA_CHR_region2 <- SHA_CHR_region[order(SHA_CHR_region$BP),]
SHA_CHR_region2[,c(4:6,8)] <- NULL

# Simple plot
ggplot(data = SHA_CHR_region2, aes(x=BP,y=LOD)) + geom_line()

# Write to file (for manual analysis)
write.table(SHA_CHR_region2, "./results/two_unit_interval/whole_region/Shannon_CHR3.txt", sep = "\t", quote = F, col.names = T, row.names = F)

max(SHA_CHR_region2$LOD) # Max LOD: 3.722

# Read file
Shannon_2USI <- read.table("./results/two_unit_interval/shannonCHR3_2USI.txt", sep = "\t", header = T)
is.unsorted(Shannon_2USI$BP)
dim(Shannon_2USI) # 177   6

# Plot the region
ggplot(Shannon_2USI, aes(x=BP,y=LOD)) + geom_line() # Looks good!

################################################################

# COMMENT

###############################################################

# So, we have 12 non-overlapping regions (Two-unit support interval)


#############################################################

#         ANALYSIS: Linkage SNPs vs Liu SNPs 

#############################################################

##### 

liu_snps_list <- readLines("./data/snpslist_liu_232.txt")

# Get all SNPs from the two-unit support intervals

All_2USI_SNPs <- list.files(path = "./results/two_unit_interval/",
                             pattern = "*.txt", full.names = TRUE)
str(All_2USI_SNPs)

# Combine all 12 interval files
linkage_snps <- do.call(rbind, lapply(All_2USI_SNPs, read.table, sep = "\t", header = T))
head(linkage_snps) 
dim(linkage_snps) # 13022     6
str(linkage_snps)

qtl_snps_list <- linkage_snps


gwas_qtl_common <- qtl_snps_list[qtl_snps_list$POS %in% liu_snps_list,]
# Write to file
write.table(gwas_qtl_common, file = "./results/SNP_overlap.txt", sep= " ", col.names = T, row.names = F, quote = F)


# Significant SNPs from pairwise similarity analysis
sig_SNP <- "rs11741861"

LODs_for_Analysis1 <- qtl_snps_list[qtl_snps_list$POS %in% sig_SNP,]

# No overlap