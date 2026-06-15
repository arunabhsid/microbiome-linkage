
#############################################################

###       GENOME WIDE QTL LINKAGE ANALYSIS

#############################################################

### MICROBIOME DATA (Quantitative Trait)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
       setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
     }

# Libraries
library(edgeR)
library(limma)
library(purrr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(data.table)
library(scales)
library(dplyr)

# Filtered out taxa with a read number of >= 10 present in 50% of samples
KIN_genera_filtered <- readRDS("./data/KIN_genera_filtered.RDS")

# Extract abundance matrix from the phyloseq object
ASV1 = as(otu_table(KIN_genera_filtered), "matrix")
# transpose if necessary
if(taxa_are_rows(KIN_genera_filtered)){ASV1 <- t(ASV1)}
# Coerce to data.frame
ASV_df = as.data.frame(ASV1)
# Get the genera names from taxonomy file
ASV_to_genus <- as.data.frame(tax_table(KIN_genera_filtered)[,6])
print(ASV_to_genus)
colnames(ASV_df) <- ASV_to_genus$Genus

# Transform these 22 taxa using log CPM
ASV_tran1 <- as.matrix(ASV_df)
# CPM transformation (from edgeR)
ASV_tran <- cpm(ASV_tran1,log = TRUE, prior.count = 1)
ASV_df_red_tr <- as.data.frame(ASV_tran)

# Visualization
# Before data transformation
jpeg("Before_transformation.jpg",height=8,width=12,units='in',res=600)
ASV_df %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()
dev.off()

# After data transformation
jpeg("After_transformation.jpg",height=8,width=12,units='in',res=600)
ASV_df_red_tr %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()
dev.off()

# Make it merge ready
microb_data <- setDT(ASV_df_red_tr, keep.rownames = TRUE)[]
colnames(microb_data)[1] <- "V2"

#############################################################

### COVARIATES

#############################################################

covariates <- readRDS("./data/covariates.RDS")

##############################################################

###         WHOLE GENOME ANALYSIS WITHOUT LD PRUNING

##############################################################

### PREPARING FINAL DATASET 

# Unlike previous analyses where we considered certain regions of
# the genome (specifically 500kb flanks on either side of the Liu risk loci),
# here we shall consider the whole genome and look for evidence of linkage

#### GENETIC DATA

#### GENERATE PED & MAP FILES USING PLINK

for (i in 1:22) {
  wg1 <- paste0("plink1 --noweb --bfile ./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort",
               " --chr ",i," --output-missing-phenotype 0 --remove ./data/indiv_remove.txt",
               " --recode --out KIN_CHR_",i)
  
  try(system(wg1))
  
}

#### PREPARE PED FILES

# Read in all the .ped files
PED_Files <- Sys.glob("KIN_CHR_*.ped")
str(PED_Files)

# Loop over all PED files
select_cols <- c("V2","V3","V4")

for (PED_file in PED_Files) {
  ped_genomewide2 <- fread(PED_file, data.table = FALSE)
  
  # Remove redundant patterns from individual ID (col 2)
  for (i in select_cols) {
    ped_genomewide2[,i] <- gsub("2019-020_pre_","",ped_genomewide2[,i])
  }
  # Merge microbiome and genetic data
  ped_micr_genomewide2 <- merge(ped_genomewide2, microb_data, 
                                by = "V2", all.x = TRUE, sort = F)[, union(names(ped_genomewide2), names(microb_data))]
  # Change NA to X
  ped_micr_genomewide2[is.na(ped_micr_genomewide2)] <- "X"
  # Merge with covariates
  ped_micr_genomewide2$AffStatus <- ped_micr_genomewide2$V6
  # FINAL MERGE
  ped_micr_cov_genomewide2 <- merge(ped_micr_genomewide2, covariates, 
                                    by = "V2", all.x = TRUE, sort = F)[, union(names(ped_micr_genomewide2), names(covariates))]
  # NAs introduced due to merge
  # Change NA to X
  ped_micr_cov_genomewide2[is.na(ped_micr_cov_genomewide2)] <- "X"
  # Write it to file (really slow!)
  write.table(ped_micr_cov_genomewide2, "KIN_CHR_temp.ped", sep = " ", quote = F, col.names = FALSE, row.names = F)
  file.rename(list.files(pattern="KIN_CHR_temp.ped"), paste0("NEW_",PED_file))
}

#### PREPARE MAP FILES

MAP_Files <- Sys.glob("KIN_CHR_*.map")
str(MAP_Files)

for(MAP_file in MAP_Files){
  
  map_genomewide2 <- fread(MAP_file, data.table = FALSE)
  map_genomewide2$V4 <- NULL
  write.table(map_genomewide2, file = "KIN_CHR_temp.map", 
              row.names = F, col.names = F, quote = F)
  file.rename(list.files(pattern = "KIN_CHR_temp.map"), paste0("NEW_",MAP_file))
  
}

#### CREATE DAT FILES

# Loop over each newly created MAP file and create DAT files 

new_map_files <- Sys.glob("NEW_KIN_CHR_*.map")
str(new_map_files)

for(new_map in new_map_files) {
  
  WG <-  paste0("awk 'BEGIN {print \"A IBD\"};{print \"M \" $2};",
                "END {print \"T Clostridium_XlVa\"};END {print \"T Roseburia\"};",
                "END {print \"T Anaerostipes\"};END {print \"T Escherichia/Shigella\"};",
                "END {print \"T Dialister\"};END {print \"T Alistipes\"};",
                "END {print \"T Bacteroides\"};END {print \"T Barnesiella\"};",
                "END {print \"T Odoribacter\"};END {print \"T Parabacteroides\"};",
                "END {print \"T Bilophila\"};END {print \"T Parasutterella\"};",
                "END {print \"T Anaerotruncus\"};END {print \"T Faecalibacterium\"};",
                "END {print \"T Subdoligranulum\"};END {print \"T Ruminococcus\"};",
                "END {print \"T Clostridium_IV\"};END {print \"T Flavonifractor\"};",
                "END {print \"T Pseudoflavonifractor\"};END {print \"T Oscillibacter\"};",
                "END {print \"T Blautia\"};END {print \"T Ruminococcus2\"};",
                "END {print \"C AffStatus\"};END {print \"C SmokStat\"};",
                "END {print \"C Age\"};END {print \"C BMI\"};\' ",
                new_map," > ",new_map,".dat")
  
  system(WG)
  
} 

### Correct extensions of these DAT files 
new_dat_files <- list.files(pattern = "*.map.dat")
str(new_dat_files)

correct_names <- gsub(".map", "", new_dat_files)
file.rename(new_dat_files, correct_names)

#### MERLIN ANALYSIS

# Run PEDSTATS
for (c in 1:22) {
  
  wg3 <- paste0("./pedstats -d NEW_KIN_CHR_",c,".dat -p NEW_KIN_CHR_",c,".ped > pedstats_KIN_CHR_",c,".txt")
  try(system(wg3))
}

# Run MERLIN REGRESS

for (d in 1:22) {
  
  ld6 <- paste0("./merlin-regress -d NEW_KIN_CHR_",d,".dat -p NEW_KIN_CHR_",d,".ped",
                " -m NEW_KIN_CHR_",d,".map --quiet --tabulate",
                " --markerNames --useCovariates --pdf --sexAsCovariate > reg_out_CHR_",d,".txt")
  try(system(ld6))
  
  # Rename pdf files to prevent overwriting
  file.rename(list.files(pattern="merlin-regress.pdf"), paste0("LOD_PLOT_CHR_",d,".pdf"))
}

### INTERPETATION GENOMEWIDE WITHOUT LD PRUNNING

AllRegressOut8 <- list.files(path = "./linkage_results_markernames_22traits/",
                             pattern = "*.tbl", full.names = TRUE)
str(AllRegressOut8)

# Combine all the 22 .tbl files
RegressFile8 <- do.call(rbind, lapply(AllRegressOut8, read.table, sep = "\t", header = T))
head(RegressFile8) 
dim(RegressFile8) # 10416802       8
str(RegressFile8)

# How many LOD scores are greater than 3?
length(which(RegressFile8$LOD > 3)) # 3481

# Sort according to LOD score
RegressResult8 <- RegressFile8[order(RegressFile8$LOD, decreasing = T)[1:3481],]

# Subset according to LOD score:
# Which quantitative traits appear the most?
print(as.data.frame(sort(table(RegressResult8$PHENOTYPE))))

# Which chromosomes appear the most (evidence of linkage)?
print(as.data.frame(sort(table(RegressResult8$CHR))))

# Look at Clostridium_XlVa
Clostridium_XlVa <- subset(RegressFile8, PHENOTYPE == "Trait: Clostridium_XlVa")
Clostridium_XlVa <- Clostridium_XlVa[order(Clostridium_XlVa$LOD, decreasing = T),]
length(which(Clostridium_XlVa$LOD > 3)) # 1644
Clostridium_XlVa[1:1644,-3]

# Look at Roseburia
Roseburia8 <- subset(RegressFile8, PHENOTYPE == "Trait: Roseburia")
Roseburia8 <- Roseburia8[order(Roseburia8$LOD, decreasing = T),]
Roseburia8[1:9,-3]
length(which(Roseburia8$LOD > 3))

# NAs

NA_SNPs_2 <- RegressFile8[which(RegressFile8$LOD == "na"),]
print(NA_SNPs_2)

# Chr 12 and 20: NAs

##############################################################

###         VISUALIZATION

##############################################################

# We need positions and not SNP name

# Run MERLIN REGRESS AGAIN
for (d in 1:22) {
  
  ld7 <- paste0("./merlin-regress -d NEW_KIN_CHR_",d,".dat -p NEW_KIN_CHR_",d,".ped",
                " -m NEW_KIN_CHR_",d,".map -t ./data/models_file.txt --quiet --tabulate",
                " --useCovariates --pdf --sexAsCovariate > reg_out_CHR_",d,".txt")
  try(system(ld7))
  
  # Rename pdf files to prevent overwriting
  file.rename(list.files(pattern="merlin-regress.pdf"), paste0("LOD_PLOT_CHR_",d,".pdf"))
}

# Results

AllRegressOut7 <- list.files(path = "./linkage_results_positions_22traits/",
                             pattern = "*.tbl", full.names = T)
str(AllRegressOut7)

# Combine all the 22 .tbl files
RegressFile7 <- do.call(rbind, lapply(AllRegressOut7, read.table, sep = "\t", header = T))
head(RegressFile7) 
dim(RegressFile7) # 10416802     8
str(RegressFile7)

# Clean up PHENOTYPE column
RegressFile7$PHENOTYPE <- gsub("Trait: ", "", RegressFile7$PHENOTYPE)

# Get the phenotypes
traits <- unique(RegressFile7$PHENOTYPE)
traits

# Checks
all_traits <- as.data.frame(table(RegressFile7$PHENOTYPE))
all_traits

# Subset data based on phenotype
trait_dfs <- as.data.frame(split(RegressFile7, RegressFile7$PHENOTYPE))
ncol(trait_dfs)

# Are POS columns identical?
identical(trait_dfs$Alistipes.POS, trait_dfs$Clostridium_XlVa.POS) # TRUE

# Plotting
as.data.frame(names(trait_dfs))
#trait_dfs_plot <- trait_dfs[,c(2,7,15,23,31,39,47,55,63,71,79,87,95)]

trait_dfs_plot <- trait_dfs[,c("Alistipes.CHR", 
                               "Alistipes.POS",
                               grep("LOD", colnames(trait_dfs), value = TRUE))]
colnames(trait_dfs_plot)[1:2] = c("CHR", "POS")

# Different chromosomes have different cM locations
as.data.frame(table(trait_dfs_plot$CHR))

######################################################################

#                     MANHATTAN PLOT

######################################################################

head(trait_dfs_plot)
str(trait_dfs_plot)

# Looking at Parasutterella
trait_dfs_plot$Parasutterella.LOD <- as.numeric(trait_dfs_plot$Parasutterella.LOD)
# Prepare a new df
trait_df_manhatt <- trait_dfs_plot[,c("CHR","POS","Parasutterella.LOD")]
head(trait_df_manhatt)

sig.dat <- trait_df_manhatt %>% 
  subset(Parasutterella.LOD > 3)

notsig.dat <- trait_df_manhatt %>% 
  subset(Parasutterella.LOD <= 3) %>%
  slice(sample(nrow(.), nrow(.) / 5))

trait_dfs_plot <- rbind(sig.dat,notsig.dat)


nCHR <- length(unique(trait_df_manhatt$CHR))
trait_df_manhatt$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(trait_df_manhatt$CHR)){
  nbp[i] <- max(trait_df_manhatt[trait_df_manhatt$CHR == i,]$POS)
  trait_df_manhatt[trait_df_manhatt$CHR == i,"BPcum"] <- trait_df_manhatt[trait_df_manhatt$CHR == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- trait_df_manhatt %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- 4
sig <- 3

jpeg("./plots/Parasutterella_manhattan.jpg",height=4,width=10,units='in',res=600)

ggplot(trait_df_manhatt, aes(x = BPcum, y = Parasutterella.LOD,
                                 color = as.factor(CHR), size = Parasutterella.LOD)) +
  geom_point(alpha = 0.75, size = 0.2) +
  geom_hline(yintercept = sig, color = "purple", linetype = "dashed", size = 0.2) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#56B4E9", "#E69F00"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome position (cM)", 
       y = "LOD Score") + 
  ggtitle("Parasutterella abundance") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, size = 6, vjust = 0.5))

dev.off()

##### VISUALIZATION: LINE PLOT TRIAL

ggplot(trait_df_manhatt, aes(x = BPcum)) +
  #geom_line(aes(y = Roseburia.LOD, color = "darkred")) +
  geom_path(aes(y = Roseburia.LOD), size = 0.2, color = "#E46726") +
  geom_hline(yintercept = sig, color = "purple", linetype = "dashed", size = 0.2) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim))  +
  labs(x = "Chromosome position (cM)", 
       y = "LOD Score") + 
  ggtitle("Regression analysis for trait: Roseburia") +
  theme_minimal() 

#############################################################

### QUANTITATIVE TRAIT LINKAGE ANALYSIS (ALPHA DIVERSITY)

#############################################################

# MICROBIOME DATA

# Load previously created phyloseq object

microbiome.genus <- readRDS("./data/KIN_genera_all.RDS")
 
# Remove genera not present in ANY samples
# Ref: https://joey711.github.io/phyloseq/plot_richness-examples.html
MG.filtered <- prune_taxa(taxa_sums(microbiome.genus) > 0, microbiome.genus)
MG.filtered

# Estimate alpha diversity
alpha.div <- estimate_richness(MG.filtered, split = T, measures = c("Chao1", "Shannon"))
alpha.div$se.chao1 <- NULL
head(alpha.div)

# Set rownames as the first column
microb.data <- data.table::setDT(alpha.div, keep.rownames = TRUE)[]
colnames(microb.data)[1] <- "V2"
str(microb.data)

# COVARIATES

covariates <- readRDS("./data/covariates.RDS")

# GENETIC DATA

# PREPARE PED FILES
# All filenames have a dot now instead of underscore to prevent overwriting

# Read in all the .ped files
PED.Files <- Sys.glob("KIN_CHR_*.ped")
str(PED.Files)

# Loop over all PED files
select.cols <- c("V2","V3","V4")

for (PED.file in PED.Files) {
  ped.genomewide2 <- data.table::fread(PED.file, data.table = FALSE)
  
  # Remove redundant patterns from individual ID (col 2)
  for (i in select.cols) {
    ped.genomewide2[,i] <- gsub("2019-020_pre_","",ped.genomewide2[,i])
  }
  # Merge microbiome and genetic data
  ped.micr.genomewide2 <- merge(ped.genomewide2, microb.data, 
                                by = "V2", all.x = TRUE, sort = F)[, union(names(ped.genomewide2), names(microb.data))]
  # Change NA to X
  ped.micr.genomewide2[is.na(ped.micr.genomewide2)] <- "X"
  # Merge with covariates
  ped.micr.genomewide2$AffStatus <- ped.micr.genomewide2$V6
  # FINAL MERGE
  ped.micr.cov.genomewide2 <- merge(ped.micr.genomewide2, covariates, 
                                    by = "V2", all.x = TRUE, sort = F)[, union(names(ped.micr.genomewide2), names(covariates))]
  # NAs introduced due to merge
  # Change NA to X
  ped.micr.cov.genomewide2[is.na(ped.micr.cov.genomewide2)] <- "X"
  # Write it to file (really slow!)
  write.table(ped.micr.cov.genomewide2, "KIN.CHR_temp.ped", sep = " ", quote = F, col.names = FALSE, row.names = F)
  file.rename(list.files(pattern="KIN.CHR_temp.ped"), paste0("AD_",PED.file))
}

#### PREPARE MAP FILES

#The previously prepared map files can be used for this analysis

#### CREATE DAT FILES

# Loop over each newly created MAP file and create DAT files 

new_map_files <- Sys.glob("NEW_KIN_CHR_*.map")
str(new_map_files)

for(new_map in new_map_files) {
  
  wgX <- paste0("awk \'BEGIN {print \"A IBD\"};{print \"M \" $2};",
                "END {print \"T Chao1\"};END {print \"T Shannon\"};",
                "END {print \"C AffStatus\"};END {print \"C SmokStat\"};",
                "END {print \"C Age\"};END {print \"C BMI\"};\' ",
                new_map," > ",new_map,".dat")
  
  system(wgX)
  
} 

### Correct extensions of these DAT files 
new_dat_files <- list.files(pattern = "*.map.dat")
str(new_dat_files)

# Careful not to overwrite previous dat files
correct_names <- gsub(".map", "", new_dat_files)
correct_names1 <- gsub("NEW","AD", correct_names)
file.rename(new_dat_files, correct_names1)

#### MERLIN ANALYSIS

# Run MERLIN REGRESS
for (d in 1:22) {
  
  ld6 <- paste0("./merlin-regress -d AD_KIN_CHR_",d,".dat -p AD_KIN_CHR_",d,".ped",
                " -m NEW_KIN_CHR_",d,".map -t --quiet --tabulate",
                " --markerNames --useCovariates --pdf --sexAsCovariate > reg_out_AlphaDiv_CHR_",d,".txt")
  try(system(ld6))
  
  # Rename pdf files to prevent overwriting
  file.rename(list.files(pattern="merlin-regress.pdf"), paste0("AD_LOD_PLOT_CHR_",d,".pdf"))
}

### INTERPETATION GENOMEWIDE (ALPHA DIVERSITY)

AllRegressOutAD <- list.files(path = ".",
                              pattern = "*.tbl", full.names = TRUE)
str(AllRegressOutAD)

# Combine all the 22 .tbl files
RegressFileAD <- do.call(rbind, lapply(AllRegressOutAD, read.table, sep = "\t", header = T))
head(RegressFileAD) 
dim(RegressFileAD) # 946982       8
str(RegressFileAD)

# How many LOD scores are greater than 3?
length(which(RegressFileAD$LOD > 3)) # 8

# Sort according to LOD score
RegressResultAD <- RegressFileAD[order(RegressFileAD$LOD, decreasing = T)[1:10],]

# Subset according to LOD score:
# Which quantitative traits appear the most?
print(as.data.frame(sort(table(RegressResultAD$PHENOTYPE))))

# Which chromosomes appear the most (evidence of linkage)?
print(as.data.frame(sort(table(RegressResultAD$CHR)))) # chromosome 3

# Look at Shannon
Shannon <- subset(RegressFileAD, PHENOTYPE == "Trait: Shannon")
Shannon <- Shannon[order(Shannon$LOD, decreasing = T),]
length(which(Shannon$LOD > 3)) # 8
Shannon_LOD3 <- Shannon[1:8,]

# MAHNATTAN PLOT 

# Clean up PHENOTYPE column
RegressFileADM <- RegressFileAD
RegressFileADM$PHENOTYPE <- gsub("Trait: ", "", RegressFileADM$PHENOTYPE)

# Get the phenotypes
traitsADM <- unique(RegressFileADM$PHENOTYPE)
traitsADM

# Checks
all_traitsADM <- as.data.frame(table(RegressFileADM$PHENOTYPE))
all_traitsADM

# Subset data based on phenotype Shannon
trait_dfsADM <- subset(RegressFileADM, PHENOTYPE == "Shannon")
ncol(trait_dfsADM)

# Plotting
trait_dfsADM_plot <- trait_dfsADM[,c("CHR", "POS", "LOD")]

# Manhattan Plot
head(trait_dfsADM_plot)
str(trait_dfsADM_plot)

# Looking at Shannon alpha diversity
sig.dat <- trait_dfsADM_plot %>% 
  subset(LOD > 3)

notsig.dat <- trait_dfsADM_plot %>% 
  subset(LOD <= 3) %>%
  slice(sample(nrow(.), nrow(.) / 5))

#trait_dfsADM_plot <- rbind(sig.dat,notsig.dat)

nCHR <- length(unique(trait_dfsADM$CHR))
trait_dfsADM_plot$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(trait_dfsADM_plot$CHR)){
  nbp[i] <- max(trait_dfsADM_plot[trait_dfsADM_plot$CHR == i,]$POS)
  trait_dfsADM_plot[trait_dfsADM_plot$CHR == i,"BPcum"] <- trait_dfsADM_plot[trait_dfsADM_plot$CHR == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- trait_dfsADM_plot %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- 4
sig <- 3

# plot
jpeg("./plots/Shannon_indexmanhattan.jpg",height=4,width=10,units='in',res=600)

ggplot(trait_dfsADM_plot, aes(x = BPcum, y = LOD,
                              color = as.factor(CHR), size = LOD)) +
  geom_point(alpha = 0.75, size = 0.2) +
  geom_hline(yintercept = sig, color = "purple", linetype = "dashed", size = 0.2) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#56B4E9", "#E69F00"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome position (cM)", 
       y = "LOD Score") + 
  ggtitle("Shannon Index") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, size = 6, vjust = 0.5))

dev.off()

#############################################################

###     EXTRACTING LINKAGE REGIONS (1-LOD-DROP)
# Superseded by the two-unit support intervals in extract_interval.R
# kept for reference, not used in the final results

#############################################################

## 1 LOD DROP from the highrest LOD score
## Arrange according to SNP-basepair position

### ROSEBURIA - CHR 19

LOD_CHR_19 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr19.tbl", 
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
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.781
# Get the first/last occurence of LOD greater than or equal to 2.781
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 2.781,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 55524813, Last occurence: BP:56578370
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 55524813, 56578370))

dim(range_of_1LODdrop) # 270   6
dim(Linkage_regions_ROS) # 332   6

# Plot the region
ggplot(data = Linkage_regions_ROS, aes(x=BP,y=LOD)) + geom_line()

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/roseburia1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### CLOSTRIDIUM_XIVa - CHR 14

LOD_CHR_14 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr14.tbl", 
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
ROS_CHR_region <- merge(LOD_Clostridium_XlVa_CHR_14, CHR_14_MAP, by = "POS", all.x = T, sort = F)
ROS_CHR_region2 <- ROS_CHR_region[order(ROS_CHR_region$BP),]
ROS_CHR_region2[,c(4:6,8)] <- NULL
# Simple plot
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 4.387
# Get the first/last occurence of LOD greater than or equal to 3.387
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 3.387,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 52758111, Last occurence: BP:88799607
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 52758111, 88799607))

dim(range_of_1LODdrop) # 864   6
dim(Linkage_regions_ROS) # 6479   6

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Clostridium_XlVa1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### CLOSTRIDIUM_XIVa - CHR 4

LOD_CHR_4 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr04.tbl", 
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
ggplot(data = ROS_CHR_region2, aes(x=BP,y=LOD)) + geom_line()
# Find the maximum LOD
max(ROS_CHR_region2$LOD) # 3.022
# Get the first/last occurence of LOD greater than or equal to 2.022
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 2.022,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 156383693, Last occurence: BP:172620583
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 156383693, 172620583))

dim(range_of_1LODdrop) # 83   6
dim(Linkage_regions_ROS) # 2475   6

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Clostridium_XlVa_CHR4_1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### Barnesiella - CHR 4

LOD_CHR_4 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr04.tbl", 
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
# Get the first/last occurence of LOD greater than or equal to 2.241
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 2.241,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 95197608, Last occurence: BP:103327536
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 95197608, 103327536))

dim(range_of_1LODdrop) # 790   6
dim(Linkage_regions_ROS) # 1062   6

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Barnesiella1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### Ruminococcus - CHR 16

LOD_CHR_16 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr16.tbl", 
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
# Get the first/last occurence of LOD greater than or equal to 3.421
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 3.421,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 58188973, Last occurence: BP:60148172
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 58188973, 60148172))

dim(range_of_1LODdrop) # 330   6
dim(Linkage_regions_ROS) # 378   6

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Ruminococcus1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### Pseudoflavonifractor - CHR 7 (CHECK IT ONCE)

LOD_CHR_7 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr07.tbl", 
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
# Get the first/last occurence of LOD greater than or equal to 2.133
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 2.133,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 35365359, Last occurence: BP:91946678
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 35365359, 91946678))

dim(range_of_1LODdrop) # 422   6
dim(Linkage_regions_ROS) # 8110   6

# Plot
ggplot(data = Linkage_regions_ROS, aes(x=BP,y=LOD)) + geom_line()

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Pseudoflavonifractor1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### Parasutterella - CHR 14

LOD_CHR_14 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr14.tbl", 
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
# Get the first/last occurence of LOD greater than or equal to 2.018
range_of_1LODdrop <- ROS_CHR_region2[ROS_CHR_region2$LOD >= 2.018,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 24792966, Last occurence: BP:32469624
# Extract all the SNPs between these base pairs
Linkage_regions_ROS <- ROS_CHR_region2 %>% 
  filter(between(BP, 24792966, 32469624))

dim(range_of_1LODdrop) # 676   6
dim(Linkage_regions_ROS) # 1234   6

# Plot
ggplot(data = Linkage_regions_ROS, aes(x=BP,y=LOD)) + geom_line()

# Write to file
write.table(Linkage_regions_ROS, "./1LOD-drop_region/Parasutterella1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### Odoribacter - CHR 22 

LOD_CHR_22 <- read.table("./linkage_results_markernames_22traits/merlin-regress-chr22.tbl", 
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
# There are two peaks. We take both separately.
# Split the dataset for the two peaks
Odori_region1 <- subset(ROS_CHR_region2, BP < 23283517)
Odori_region2 <- subset(ROS_CHR_region2, BP > 23283517)

# Find the maximum LOD for both the regions
max(Odori_region1$LOD) # 4.311
max(Odori_region2$LOD) # 4.259

# REGION 1
# Get the first/last occurence of LOD greater than or equal to 3.311
range_of_1LODdrop_1 <- Odori_region1[Odori_region1$LOD >= 3.311,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop_1$BP)
# First occurence: BP: 22782258, Last occurence: BP:22960284
# Extract all the SNPs between these base pairs
Linkage_regions_Odori_1 <- Odori_region1 %>% 
  filter(between(BP, 22782258, 22960284))

dim(range_of_1LODdrop_1) # 10   6
dim(Linkage_regions_Odori_1) # 10   6

# Write to file
write.table(Linkage_regions_Odori_1, "./1LOD-drop_region/Odoribacter1LODdrop_region1.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# REGION 2
# Get the first/last occurence of LOD greater than or equal to 3.259
range_of_1LODdrop_2 <- Odori_region2[Odori_region2$LOD >= 3.259,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop_2$BP)
# First occurence: BP: 31487354, Last occurence: BP:37013694
# Extract all the SNPs between these base pairs
Linkage_regions_Odori_2 <- Odori_region2 %>% 
  filter(between(BP, 31487354, 37013694))

dim(range_of_1LODdrop_2) # 572   6
dim(Linkage_regions_Odori_2) # 1139   6

# Write to file
write.table(Linkage_regions_Odori_2, "./1LOD-drop_region/Odoribacter1LODdrop_region2.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###

### ALPHA DIVERSITY SHANNON INDEX - CHR 3

LOD_CHR_3 <- read.table("./linkage_results_AD_markernames/merlin-regress-chr03.tbl", 
                        sep = "\t", header = T)
LOD_Shannon_CHR_3 <- subset(LOD_CHR_3, PHENOTYPE == "Trait: Shannon")
dim(LOD_Shannon_CHR_3) # 32057     8
LOD_Shannon_CHR_3$PHENOTYPE <- gsub("Trait:", "", LOD_Shannon_CHR_3$PHENOTYPE)
max(LOD_Shannon_CHR_3$LOD) # Max LOD: 3.722

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
# Get the first/last occurence of LOD greater than or equal to 2.722
range_of_1LODdrop <- SHA_CHR_region2[SHA_CHR_region2$LOD >= 2.722,]
# Check if the region is not sorted
is.unsorted(range_of_1LODdrop$BP)
# First occurence: BP: 189910072, Last occurence: BP:189981934
# Extract all the SNPs between these base pairs
Linkage_regions_SHA <- SHA_CHR_region2 %>% 
  filter(between(BP, 189910072, 189981934))

dim(range_of_1LODdrop) # 13  6
dim(Linkage_regions_SHA) #2131   6

#Plot
ggplot(data = Linkage_regions_SHA, aes(x=BP,y=LOD)) + geom_line()
# Write to file
write.table(Linkage_regions_SHA, "1LOD-drop_region/shannon1LODdrop.txt", sep = "\t", quote = F, col.names = T, row.names = F)




