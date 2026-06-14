#############################################################

### HERITABILITY ANALYSIS (AT GENUS LEVEL)

#############################################################
# Run this script from the folder that contains it (the folder with the ./data subfolder).

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Libraries
library(edgeR)
library(limma)
library(purrr)
library(tidyr)
library(ggplot2)
library(phyloseq)

### Microbiome data preparation

microb.genus.filtered <- readRDS("./data/KIN_genera_filtered.RDS")
# Filtered out taxa with a read number of >= 10 present in 50% of samples
# OR, Removed taxa not seen more than 10 times in at least 50% of the samples

# Extract ASV table from phyloseq obj
ASV.table = as(otu_table(microb.genus.filtered), "matrix")
# transpose if necessary
if(taxa_are_rows(microb.genus.filtered)){ASV.table <- t(ASV.table)}
# Coerce to data.frame
ASV.table = as.data.frame(ASV.table)
# Get the genera names from taxonomy file
ASV.to.genus <- as.data.frame(tax_table(microb.genus.filtered)[,6])
print(ASV.to.genus)
colnames(ASV.table) <- ASV.to.genus$Genus

# Transform these 22 taxa using log CPM
genera.transf <- as.matrix(ASV.table)
# CPM transformation (from edgeR)
genera.transformed <- cpm(genera.transf,log = TRUE, prior.count = 1)
genera.transformed <- as.data.frame(genera.transformed)

# Plot
genera.transformed %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# Make rownames as the first column (for merging)
gen.tran.final <- data.table::setDT(genera.transformed, keep.rownames = T)[]
colnames(gen.tran.final)[1] <- "V2"

########################################################

###              GENOTYPE DATA
# Generate .ped and .map files for a dummy SNP

########################################################

# Generate ped and map files using Plink
liu <- paste0("plink1 --noweb --bfile ./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort",
             " --output-missing-phenotype 0 --remove ./data/indiv_remove.txt",
             " --extract ./data/dummy_snp.txt --recode --out Heritability_file")

system(liu)

# Read in .ped file
ped_liuwide <- data.table::fread("Heritability_file.ped", data.table = F)
dim(ped_liuwide)  # 1211 8
head(ped_liuwide[1:5,1:8])

# Remove redundant patterns 
select_cols <- c("V2","V3","V4")
for (i in select_cols) {
  ped_liuwide[,i] <- gsub("2019-020_pre_","",ped_liuwide[,i])
}
head(ped_liuwide[1:5,1:8])

# Clean up .map files
map_liuwide <- read.table("Heritability_file.map", sep = "")
nrow(map_liuwide) # 1 SNP
map_liuwide$V4 <- NULL
write.table(map_liuwide, file = "Heritability_file.map", 
            row.names = F, col.names = F, quote = F)

#### MERGE WITH MICROBIOME DATA

dim(ped_liuwide) #1211 8
dim(gen.tran.final) #940 33
head(ped_liuwide)[1:5,1:8]
head(gen.tran.final)

# Merge microbiome and genetic data
ped_micr_liuwide <- merge(ped_liuwide, gen.tran.final, by = "V2", all.x = TRUE, sort = F)
dim(ped_micr_liuwide) #1211 30
# NAs introduced due to merge
# Change NA to X
ped_micr_liuwide[is.na(ped_micr_liuwide)] <- "X" 

head(ped_micr_liuwide[1:5,1:10])
# Reorder V1 to first position
ncol(ped_micr_liuwide) # 30
ped_micr_liuwide <- ped_micr_liuwide[,c(2,1,3:30)]

#### PREPARE COVARIATES

# Other Covariates (excluding Household)
cov.heritability <- readRDS("./data/covariates.RDS")
head(cov.heritability)

#### MERGE WITH COVARIATES

# Affection Status (V6 is the affection status)
ped_micr_liuwide$AffStatus <- ped_micr_liuwide$V6

# FINAL MERGE
ped_micr_cov_liuwide <- merge(ped_micr_liuwide, cov.heritability, by = "V2", all.x = TRUE, sort = F)
dim(ped_micr_cov_liuwide) #1211 34
# NAs introduced due to merge
# Change NA to X
ped_micr_cov_liuwide[is.na(ped_micr_cov_liuwide)] <- "X"
head(ped_micr_cov_liuwide[1:5,1:10])
# Reoder V1 to first position
ncol(ped_micr_cov_liuwide) #34
ped_micr_cov_liuwide <- ped_micr_cov_liuwide[,c(2,1,3:34)]

# Write it to file (really slow!)
write.table(ped_micr_cov_liuwide, "ped_kin_liuwide_heritability.ped", sep = " ", quote = F, col.names = FALSE, row.names = F)

### Create the lines of code for .dat file using awk
# Names of taxa to be tested for heritability
head(ASV.to.genus)
iterator <- as.character(ASV.to.genus$Genus)

#for (genus in iterator) {
  #code.line2 <- paste0("END {print \"T ",iterator,"\"};", collapse = "")

#}
code.line2 <- paste0("END {print \"T ",iterator,"\"};", collapse = "")
code.line1 <- paste0("awk \'BEGIN {print \"A IBD\"};{print \"M \" $2};")
code.line3 <- paste0("END {print \"C AffStatus\"};END {print \"C SmokStat\"};",
                     "END {print \"C Age\"};END {print \"C BMI\"};\'",
                     " Heritability_file.map > Heritability_file.dat")
code.lines <- paste0(code.line1,code.line2,code.line3)
code.lines

#### CREATE DAT FILE

dat.call <- code.lines
try(system(dat.call))

#### Run MERLIN for Variance Components

merlin.call <- paste0("./merlin -d Heritability_file.dat",
                      " -p ped_kin_liuwide_heritability.ped",
                      " -m Heritability_file.map --vc --singlepoint --quiet",
                      " --useCovariates --sexAsCovariate > VC_her_22_traits_final_v2.txt")

try(system(merlin.call))


# Results

herit_result <- readLines("VC_her_22_traits_final_v2.txt")
h2 = grep("Phenotype", herit_result)
writeLines(herit_result[h2], con = paste0("heritability_result_22_traits.txt"))

# Read the result file
herit_result2 <- read.table("heritability_result_22_traits.txt", sep = "")
herit_result3 <- dplyr::select(herit_result2, V2, V8)
data.table::setnames(herit_result3, old = c("V2","V8"), new = c("Genus","h2"))
herit_result3$h2 <- gsub("%)","",herit_result3$h2)
herit_result3$h2 <- as.numeric(herit_result3$h2)
herit_result3$Genus <- as.factor(herit_result3$Genus)
herit_res_final <- herit_result3[order(herit_result3$h2),]

# Table for paper
write.table(herit_res_final,"heritability_table.txt", sep = ",", quote = F,row.names = F, col.names = T)

######################## Final plots ################### 

# order the factor levels by h2 values
str(herit_result3)
herit_result3$Genus <- factor(herit_result3$Genus, levels = herit_result3$Genus[order(herit_result3$h2)])
herit_result3$Genus

pdf("Heritability.pdf")
ggplot(herit_result3, aes(x = h2, y = Genus)) + theme_minimal() + geom_bar(stat = "identity") +
  labs(x = "h² (%)")
dev.off()

tiff("Heritability_tiff.tiff",height=4,width=6,units='in',res=600)
ggplot(herit_result3, aes(x = h2, y = Genus)) + theme_minimal() + geom_bar(stat = "identity") +
  labs(x = "h² (%)")
dev.off()

jpeg("Heritability_jpeg.jpg",height=4,width=6,units='in',res=600)
ggplot(herit_result3, aes(x = h2, y = Genus)) + theme_minimal() + geom_bar(stat = "identity") +
  labs(x = "h² (%)")
dev.off()
