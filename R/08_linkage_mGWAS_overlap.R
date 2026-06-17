#################################################################

# OVERLAP OF LINKAGE REGIONS WITH PREVIOUSLY REPORTED mGWAS SNPS 

#################################################################


if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
       setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
     }

# Libraries
library(dplyr)
library(tidyr)

# Get the two-unit support interval files

All_intervals_SNPs <- list.files(path = "./results/two_unit_interval/",
                                 pattern = "*.txt", full.names = TRUE)
str(All_intervals_SNPs)

# Combine all the 22 .tbl files
linkage_snps <- do.call(rbind, lapply(All_intervals_SNPs, read.table, sep = "\t", header = T))
head(linkage_snps) 
dim(linkage_snps) # 13022     6
str(linkage_snps)

# Change names
qtl_snps_list <- linkage_snps

mGWAS_snps_list <- readLines("./data/mGWAS_SNPs_all.txt")
mgwas_qtl_common <- qtl_snps_list[qtl_snps_list$POS %in% mGWAS_snps_list,] # 3 common SNPs found

# Are all mGWAS SNPs present in our data?

All_SNPS_Kindred <- read.table("./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort.bim")
All_SNPS_Kindred$V2 <- gsub("GSA-","", All_SNPS_Kindred$V2)

length(mGWAS_snps_list) # 379

common_snps <- All_SNPS_Kindred[All_SNPS_Kindred$V2 %in% mGWAS_snps_list,]
common_snps1 <- common_snps$V2
ungenotyped_snps <- setdiff(mGWAS_snps_list, common_snps1) # 70

## BUBIER SNPS (2021- review and contains all previous SNPs) - That are present in chromosomes 3,4,7,14,16,19,22
bubier_snps <- read.table("./data/bubier_snps.csv", sep =",", header = T)
dim(bubier_snps) # 52   2
length(unique(bubier_snps$SNP)) # 49
bubs <- bubier_snps$SNP
kins <- All_SNPS_Kindred$V2
# Are Bubier SNPS present in  Kindred
common_snps_bub <- All_SNPS_Kindred[All_SNPS_Kindred$V2 %in% bubier_snps$SNP,] # Only 4

# Find the positions in VEP as many of Bubier SNPs were not genotyped

bub_positions <- read.table("./data/bubier_snps_position.txt")
bub_pos <- bub_positions[,1:2]
bub_pos1 <- bub_pos[!duplicated(bub_pos),]
bub_pos2 <- bub_pos1 %>% tidyr::separate(V2, into = c("CHR","POS"), sep = ":")
bub_pos3 <- bub_pos2 %>% separate(POS, into = c("pos_from","pos_to"), sep = "-")
bub_pos3$pos_from <- as.numeric(bub_pos3$pos_from)
bub_pos3$pos_MB <- bub_pos3$pos_from/1000000

# Use bub_pos3 to find the overlap

# bubier/mGWAS SNPs are largely ungenotyped here, so the rs-ID match returns nothing.
# match by genomic position instead: for each interval, keep the bubier SNPs whose
# position falls inside the interval's bp range.
bub_pos3$CHR <- as.numeric(bub_pos3$CHR)

bubier_overlap <- do.call(rbind, lapply(All_intervals_SNPs, function(f) {
  iv <- read.table(f, sep = "\t", header = TRUE)
  hits <- bub_pos3[bub_pos3$CHR == iv$CHR[1] &
                     bub_pos3$pos_from >= min(iv$BP) &
                     bub_pos3$pos_from <= max(iv$BP), ]
  if (nrow(hits) > 0) { hits$interval <- basename(f); hits } else NULL
}))

bubier_overlap
write.table(bubier_overlap, "./results/bubier_overlap.txt", sep = "\t", quote = F, row.names = F)