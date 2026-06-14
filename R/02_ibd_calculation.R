
#############################################################

### IDENTITY BY DESCENT (IBD) CALCULATION

#############################################################

# Path
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
       setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
     }

# Libraries
library(foreach)
library(ggplot2)
library(data.table)
library(pedtools)
library(ribd)
library(kinship2)
library(tidyr)

# List of Liu SNPs
input_names <- readLines("./data/snpslist_liu_232.txt")
str(input_names)
head(input_names)

########################################################

### Generate .ped and .map files for each of the Liu SNPs 
### with a window size of 500kb 

########################################################

foreach(i = input_names) %do% {
  
  cmd <-  paste0("plink1 --noweb --bfile ./data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort",
  " --output-missing-phenotype 0 --remove ./data/indiv_remove.txt --snp ", i, " --window 500",
  " --recode --out IBD_calc/ibd_", i)
  system(cmd)
  
}

# QC: number of SNPs in each map file (500 kb region per SNP)
# Total SNPs: 204

# Read in all the .map files
map_Files <- Sys.glob("IBD_calc/*.map")
str(map_Files)

# Loop over all the map files
for (map_file in map_Files) {
  
  # read map file
  map <- read.table(map_file, sep = "")
  
  # create a new dataframe
  SNP <- data.frame(File = map_file, SNP_num = nrow(map))
  
  # write new data to separate file:
  write.table(SNP, 
              "IBD_calc/results/SNP_density.csv",
              append = TRUE,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
}

# Read in the file
f1 <- read.csv("IBD_calc/results/SNP_density.csv", header = FALSE)
head(f1)

# Clean up SNP names
f1$V1 <- gsub("IBD_calc/ibd_|.map","", f1$V1)
head(f1)
min(f1$V2) # 65 SNPs
max(f1$V2) # 493 SNPs

# Create a plot
jpeg("IBD_calc/results/Number of SNPs at window size 500kb.jpg",height=8,width=12,units='in',res=600)
ggplot(f1, aes(x = V1, y = V2)) + geom_col() + coord_flip() + 
  labs(title = "Number of SNPs at window size 500kb", x = "Map file", y = "Number of SNPs") +
  theme(axis.text.y = element_text(size = 4)) 
dev.off()

# Histogram
jpeg("IBD_calc/results/Number of SNPs at window size 500kb_hist.jpg",height=8,width=12,units='in',res=600)
p<-ggplot(f1, aes(x=V2)) + 
  geom_histogram(color="black", fill="grey") + 
  labs(title = "Number of SNPs at window size 500kb", x = "Map file of size", y = "Number of SNPs")
p+ geom_vline(aes(xintercept=mean(V2)),
              color="blue", linetype="dashed", size=1)
dev.off()

########################################################

### Generate .dat files for every .map file (using awk)

########################################################

# Mapfiles
mapfiles <-Sys.glob(path = "IBD_calc/*.map")
str(mapfiles)
# Create dat files
foreach(i = mapfiles) %do% {
  cmd <- paste0("awk \'BEGIN {print \"A IBD\"};{print \"M \" $2}\' ",i," > ",i,".dat")
  system(cmd)
}

########################################################

### Renaming .dat files

########################################################

datfiles <- Sys.glob(path = "IBD_calc/*.dat")
head(datfiles)

new_names <- gsub(".map", "", datfiles)
file.rename(datfiles, new_names)

########################################################

### Cleaning up the .map files

########################################################

# Read in all the .map files
map_Files2 <- Sys.glob(path = "IBD_calc/*.map")
str(map_Files2)

# Loop over all the map files
for (map_file in map_Files2) {
  
  # read map file
  map <- read.table(map_file, sep = "")
  head(map)
  # remove the fourth column (basepairs)
  map$V4 <- NULL
  
  # Overwrite each file
  write.table(map, 
              map_file,
              append = FALSE,
              quote = FALSE,
              sep = " ",
              row.names = FALSE,
              col.names = FALSE)
}

########################################################

### IBD calculation with MERLIN

########################################################

setwd("IBD_calc/")
# We could us the SNP names as iterator
head(input_names)

# Merlin can also be installed globally using make file
foreach(i = input_names) %do% {
  cmd <- paste0("./merlin -d ibd_",i,".dat -p ibd_",i,".ped -m ibd_",i,".map --markerNames --ibd > MER_IBD_",i,".txt")
  system(cmd)
  file.rename(list.files(pattern="merlin.ibd"), paste0("IBD_",i,".ibd"))
}

########################################################

### Check the .ibd files

########################################################

# Read in all the .map files
ibd_Files <- Sys.glob("*.ibd")
str(ibd_Files)

# Loop over all the ibd files
for (ibd_file in ibd_Files) {
  
  # read ibd file
  ibd <- read.table(ibd_file, sep = "")
  
  # create new dataframe
  Row <- data.frame(File = ibd_file, row_num = nrow(ibd))
  
  # write new data to separate file:
  write.table(Row, 
              "./results/IBD_result_length.csv",
              append = TRUE,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
}

# Read in the file
f2 <- read.csv("./results/IBD_result_length.csv", header = FALSE)
head(f2)

# Clean up SNP names
f2$V1 <- gsub("IBD_|.ibd","", f2$V1)
head(f2)
min(f2$V2) # 159614
max(f2$V2)
# Plot
library(ggplot2) # 1588345

# Plot
jpeg("./results/Number of rows per ibd file.jpg",height=8,width=12,units='in',res=600)
ggplot(f2, aes(x = V1, y = V2)) + geom_col() + coord_flip() + 
  labs(title = "Number of rows per .ibd file", x = "IBD file", y = "Number of Rows") +
  theme(axis.text.y = element_text(size = 4)) 
dev.off()

# QC: Check the max of P0,P1 and P2 for each pair 
# and check that we have more values close to 1 in the histogram

# Read all ibd files
ibd_Files <- Sys.glob("*.ibd")
str(ibd_Files)

# Loop over all the ibd files
for (ibd_file in ibd_Files) {
  
  ibd <- read.table(ibd_file, sep = "")
  
  # subset each ibd file keeping only P0,P1,P2
  ibd1 <- ibd[,5:7]
  
  # calculate row max
  max <- apply(X=ibd1, MARGIN=1, FUN=max)
  
  # create new dataframe
  My_df <- data.frame(row_max = max)
  
  # write new data to separate file
  write.table(My_df, 
              "./results/Row_max.csv",
              append = TRUE,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
}

# Read in the file
f3 <- read.csv("./results/Row_max.csv", header = FALSE)
nrow(f3) # 145969963
min(f3$V1) # 0.33333
max(f3$V1) # 1

# Plot
jpeg("./results/Max_IBD_all_pairs.jpg",height=8,width=12,units='in',res=600)
f3$V1 <- as.numeric(f3$V1)
hist(f3$V1)
dev.off()

# Keep only Liu SNPs

# Cleaning up each of the 204 .ibd files before rbinding
# Select only rows containing Liu SNPs (discard the flanking SNPs)

# Read in all the .ibd files
IBDFILES <- list.files(pattern = "*.ibd")
str(IBDFILES)

# Loop over all the ibd files
for (IBDFILE in IBDFILES) {
  
  # read ibd files
  IBDF <- read.table(IBDFILE, sep = "")
  head(IBDF)
  
  # KEEP ONLY LIU SNPs
  
  # clean up the filename
  cleaner <- gsub("IBD_|.ibd","",IBDFILE)
  # keep rows which have cleaner in them
  IBDF2 <- IBDF[IBDF$MARKER %in% cleaner,]
  
  # Overwrite each file
  write.table(IBDF2, 
              IBDFILE,
              append = FALSE,
              quote = FALSE,
              sep = " ",
              row.names = FALSE,
              col.names = FALSE)
}

cat("---------- IBD calculation complete ----------\n")

########################################################

### Filtering out uninformative pairs

########################################################

### if re-reun needed, you can start directly from here!

# Make a long table combining all the 204 .ibd tables
ibd_out_files <- list.files(pattern = "*.ibd")

allIBDFiles <- do.call(rbind, lapply(ibd_out_files, fread, data.table = F, header = FALSE))
head(allIBDFiles)
dim(allIBDFiles)   #735301      7

# Set new column names
setnames(allIBDFiles, old = c("V1","V2","V3","V4","V5","V6","V7"),
         new = c("FAMILY","ID1","ID2","MARKER","P0","P1","P2"))

# Remove redundant pattern from genotyped IDs
sapply(allIBDFiles, class) # check class of each var
allIBDFiles$ID1 <- gsub("2019-020_pre_","",allIBDFiles$ID1)
allIBDFiles$ID2 <- gsub("2019-020_pre_","",allIBDFiles$ID2)

# Now, we can proceed with the filtering
head(allIBDFiles)

# Delete self pairs
allIBDFiles2 <- subset(allIBDFiles, ID1 != ID2)
nrow(allIBDFiles2) # 506537
head(allIBDFiles2)

# Delete rows when there is a un-genotyped ID (*_*) 
## Individuals with a 6 character tokens are those, which where previously the 1_2 etc. 
allIBDFiles3 <- allIBDFiles2[!grepl("\\b[a-zA-Z0-9]{6}\\b", allIBDFiles2$ID1),]
allIBDFiles4 <- allIBDFiles3[!grepl("\\b[a-zA-Z0-9]{6}\\b", allIBDFiles3$ID2),]
nrow(allIBDFiles4) # 198587
head(allIBDFiles4)

# Parent-Child + MZ + Unrelated pairs filtering

# FIRST CALCULATE THE KAPPAS
# Read the .fam file
fam_file <- read.table("../data/KIN_GSA_MD_24v2_0_B1_b37_qc_completed_sort.fam")
head(fam_file)
# Remove sex and phenotype column
fam_file$V5 <- NULL
fam_file$V6 <- NULL
head(fam_file)
# Set column names
colnames(fam_file) <- c("pedigree_id", "id","fatherID", "motherID")
# Remove redundant patterns
fam_file$id <- gsub("2019-020_pre_","",fam_file$id)
fam_file$fatherID <- gsub("2019-020_pre_","",fam_file$fatherID)
fam_file$motherID <- gsub("2019-020_pre_", "", fam_file$motherID)
head(fam_file)

# Write to a file (It's already written)
write.table(fam_file, "new_pedigree.txt", sep = " ", row.names = FALSE, quote = F)

# Read all pedigrees into a list
x = readPed("new_pedigree.txt", header = T, famid = 1, id = 2, fid = 3, mid = 4)

# Compute kappa coefficients of all pairs in each pedigree
kappas = lapply(x, kappaIBD)

# All Parent-Child + MZ + Unrelated pairs in each pedigree
PO = lapply(kappas, function(df) subset(df, kappa1 == 1, select = c(id1,id2)))
MZ = lapply(kappas, function(df) subset(df, kappa0 == 0 & kappa1 == 0 & kappa2 == 1, select = c(id1,id2)))
UN = lapply(kappas, function(df) subset(df, kappa0 == 1 & kappa1 == 0 & kappa2 == 0, select = c(id1,id2)))
# All sib-pairs pairs in each pedigree
SS = lapply(kappas, function(df) subset(df, kappa0 == 0.25 & kappa1 == 0.50 & kappa2 == 0.25, select = c(id1,id2)))
# Half sibling, Grandfather-grandchild, Uncle-Nephew
HS = lapply(kappas, function(df) subset(df, kappa0 == 0.5 & kappa1 == 0.50 & kappa2 == 0, select = c(id1,id2)))
# Double first cousins
DFC = lapply(kappas, function(df) subset(df, kappa0 == 0.5625 & kappa1 == 0.375 & kappa2 == 0.0625, select = c(id1,id2)))
# First cousins
FC = lapply(kappas, function(df) subset(df, kappa0 == 0.75 & kappa1 == 0.25 & kappa2 == 0, select = c(id1,id2)))

# Combine them into a single dataframe
PO_all = do.call(rbind, PO)
# MZ_all = do.call(rbind, MZ) # Can't find MZ pairs wihout genotype info.
UN_all = do.call(rbind, UN)
SS_all = do.call(rbind, SS)
HS_all = do.call(rbind, HS)
DFC_all = do.call(rbind, DFC) # No double first cousins
FC_all = do.call(rbind, FC)

head(PO_all)
head(UN_all)
head(SS_all)

# Create pairs: concatenate id1 and id2 from SS_all
SS_all$PAIR <- ifelse(SS_all$id1<SS_all$id2, paste0(SS_all$id1,"-",SS_all$id2), paste0(SS_all$id2,"-",SS_all$id1))
SS_all$Sib_pair <- SS_all$PAIR
head(SS_all)
# Concatenate id1 and id2 from PO_all
PO_all$PAIR <- ifelse(PO_all$id1<PO_all$id2, paste0(PO_all$id1,"-",PO_all$id2), paste0(PO_all$id2,"-",PO_all$id1))
head(PO_all)
# Concatenate id1 and id2 from UN_all
UN_all$PAIR <- ifelse(UN_all$id1<UN_all$id2, paste0(UN_all$id1,"-",UN_all$id2), paste0(UN_all$id2,"-",UN_all$id1))
head(UN_all)
# Similarly, do this for all
HS_all$PAIR <- ifelse(HS_all$id1<HS_all$id2, paste0(HS_all$id1,"-",HS_all$id2), paste0(HS_all$id2,"-",HS_all$id1))
head(HS_all)
FC_all$PAIR <- ifelse(FC_all$id1<FC_all$id2, paste0(FC_all$id1,"-",FC_all$id2), paste0(FC_all$id2,"-",FC_all$id1))
head(FC_all)
# MZ Twins (refer mz_twins.txt)
MZ_all <- data.frame("PAIR" = c("Ao9aa-Moo8o","Moo8o-Ao9aa","FoZo7-Ii2fa","Ii2fa-FoZo7",
                                "aeY4o-caiH5","caiH5-aeY4o","auY1i-uZie3","uZie3-auY1i",
                                "eiCh8-aer7K","aer7K-eiCh8","mai7I-tee7E","tee7E-mai7I"))
MZ_all$PAIR <- as.character(MZ_all$PAIR)
head(MZ_all)

########################################################

###                      FILTERING

########################################################

dim(allIBDFiles4)  # 198587       7
dim(PO_all)   #1238    2

# Concatenate ID1 and ID2 from allIBDFile4
allIBDFiles4$PAIR <- ifelse(allIBDFiles4$ID1<allIBDFiles4$ID2, paste0(allIBDFiles4$ID1,"-",allIBDFiles4$ID2), paste0(allIBDFiles4$ID2,"-",allIBDFiles4$ID1))
head(allIBDFiles4)

# Finally, remove all pairs contained in PO_all from allIBDFile4
allIBDFiles5 <- allIBDFiles4[!allIBDFiles4$PAIR %in% PO_all$PAIR,]
dim(allIBDFiles5)   #119485      9
head(allIBDFiles5)

# Finally, remove all pairs contained in UN_all from allIBDFile5
allIBDFiles6 <- allIBDFiles5[!allIBDFiles5$PAIR %in% UN_all$PAIR,]
dim(allIBDFiles6)   #99096      9
head(allIBDFiles6)

# Finally, remove all pairs contained in MZ_all from allIBDFile6
allIBDFiles7 <- allIBDFiles6[!allIBDFiles6$PAIR %in% MZ_all$PAIR,]
dim(allIBDFiles7)   #97962     9
head(allIBDFiles7)

########################################################

# QC: Check the distribution of max values
ibd_result_now <- allIBDFiles7[,5:7]

ibd_result_now$max<-apply(X=ibd_result_now, MARGIN=1, FUN=max)
str(ibd_result_now)

ibd_result_now$max <- as.numeric(ibd_result_now$max)
hist(ibd_result_now$max)
### Looks good

# Calculate expected value and append it to dataframe
columns <- names(allIBDFiles7)[5:7]
allIBDFiles7[columns] <- lapply(allIBDFiles7[columns], as.numeric)
allIBDFiles7$EXP.VAL <- (0*allIBDFiles7$P0+1*allIBDFiles7$P1+2*allIBDFiles7$P2)
head(allIBDFiles7)
dim(allIBDFiles7) # 97962     9
# So, allIBDFiles7 has all the pairs except UN, PO and MZ pairs

########################################################

###       UNWEIGHTED SCORE (JUMP TO: ANALYSIS PER SNP)

########################################################

# ADD UP the Exp. value for all the Liu SNPs for each pair
allIBDFiles8 <- allIBDFiles7[,c(8,9)]
head(allIBDFiles8)

allIBDFiles9 <- aggregate(EXP.VAL ~ PAIR, data = allIBDFiles8, FUN = sum)
head(allIBDFiles9)
dim(allIBDFiles9)  # 519   2

# QC: Checking relationship pairs

# relationships from the pedigree file
dim(SS_all) # Sibling pairs: 384   4
dim(FC_all) # First cousins: 16  3
dim(HS_all) # Half sibling + GG + AV: 392   3

# Our data
dim(allIBDFiles9)   #519   2

# So, 519 rows should be: SS+FC+HS. Let's saggregate them.
sibling_pairs <- allIBDFiles9[allIBDFiles9$PAIR %in% SS_all$PAIR,]
dim(sibling_pairs)   #352   2
first_cousin_pairs <- allIBDFiles9[allIBDFiles9$PAIR %in% FC_all$PAIR,]
dim(first_cousin_pairs)   #4   2
half_sibling_pairs <- allIBDFiles9[allIBDFiles9$PAIR %in% HS_all$PAIR,]
dim(half_sibling_pairs)   #158   2

# sibling_pairs,first_cousin_pairs,half_sibling_pairs add up to 514,
# which means that there are 5 pairs which are of other relationship type
# I checked them. They belong to Family 9 (Old: Family 246)
# Their expected probabilities can be found from kappaIBD values 
# of the ped file (WITHOUT GENOTYPES). Just use them.

all_pairs <- rbind(sibling_pairs,first_cousin_pairs,half_sibling_pairs)
missing_ids <- allIBDFiles9[!allIBDFiles9$PAIR %in% all_pairs$PAIR,]

# The missing ids belong to family 9

###########################################

# Plotting family 9
fam9 <- as.data.frame(x["9"])
fam9$ped <- "9"
head(fam9)
fam9 <- fam9[c(5,1,2,3,4)]
setnames(fam9, old = c("ped","X9.id","X9.fid","X9.mid","X9.sex"),
         new = c("ped","id","father","mother","sex"))
head(fam9)
str(fam9)
fam9$sex[fam9$sex == 0] <- 2 # Zeroes don't work
fam9$father[fam9$father == 0] <- "" # NA, 0 doesn't wok. Weird. So, empty cell.
fam9$mother[fam9$mother == 0] <- ""
# Ped object (list)
pedAll <- pedigree(id = fam9$id, dadid = fam9$father, momid = fam9$mother,
                   sex = fam9$sex, famid = fam9$ped)
print(pedAll)
# Plot family 9
jpeg("./results/Family_9.jpg",height=8,width=12,units='in',res=600)
plot(pedAll["9"])
dev.off()

# Get the estimated values from family 9 (kappas df)
family_9 <- as.data.frame(kappas$`9`)
family_9$PAIR <- ifelse(family_9$id1<family_9$id2, paste0(family_9$id1,"-",family_9$id2), paste0(family_9$id2,"-",family_9$id1))
family_9$E_exp <- (0*family_9$kappa0+1*family_9$kappa1+2*family_9$kappa2)
family_9_final <- family_9[,c(6,7)]  
# merge
other_pairs <- merge(missing_ids, family_9_final, by = "PAIR")

########################################################

###              ANALYSIS PER SNP

########################################################
head(allIBDFiles7)

allIBDFiles100 <- allIBDFiles7[,c(8,4,9)]
dim(allIBDFiles100)   # 97962     3
setnames(allIBDFiles100, old = "EXP.VAL", new = "E_obs")
head(allIBDFiles100)

########################################################

###           STANDARDIZED IBD VALUES 

########################################################

# CREATING A COVARIATE "std_ibd"
# Scaling (groupwise mean subtraction)
# Segregate different types of relationships (pairs)
# So, 519 rows should be: SS+FC+HS. Let's segregate them.

sibling_pairs_2 <- allIBDFiles100[allIBDFiles100$PAIR %in% SS_all$PAIR,]
dim(sibling_pairs_2)   # 66444     3
first_cousin_pairs_2 <- allIBDFiles100[allIBDFiles100$PAIR %in% FC_all$PAIR,]
dim(first_cousin_pairs_2)   # 756   3
half_sibling_pairs_2 <- allIBDFiles100[allIBDFiles100$PAIR %in% HS_all$PAIR,]
dim(half_sibling_pairs_2)   # 29817     3

# Find the missing pairs
all_pairs2 <- rbind(sibling_pairs_2,first_cousin_pairs_2,half_sibling_pairs_2) #97017 rows
missing_pairs <- allIBDFiles100[!allIBDFiles100$PAIR %in% all_pairs2$PAIR,] # 945 rows
# So, now it adds up to 97962 pairs!

# Add E_exp column to the 3 dfs + remaining 5 pairs
sibling_pairs_2$E_exp <- 1
first_cousin_pairs_2$E_exp <- 0.25
half_sibling_pairs_2$E_exp <- 0.5

# Family_9_final has the expected ibd values already
other_pairs_2 <- merge(missing_pairs, family_9_final, by = "PAIR")

# Add "std_ibd" variable (E_obs - E_exp)
sibling_pairs_2$std_ibd <- sibling_pairs_2$E_obs - sibling_pairs_2$E_exp
first_cousin_pairs_2$std_ibd <- first_cousin_pairs_2$E_obs - first_cousin_pairs_2$E_exp
half_sibling_pairs_2$std_ibd <- half_sibling_pairs_2$E_obs - half_sibling_pairs_2$E_exp
other_pairs_2$std_ibd <- other_pairs_2$E_obs - other_pairs_2$E_exp

# Join all the four dfs
std_ibd_df <- rbind(sibling_pairs_2, first_cousin_pairs_2, half_sibling_pairs_2, other_pairs_2)
nrow(std_ibd_df)   # 97962
head(std_ibd_df)
std_ibd_df_final <- std_ibd_df[,-c(3,4)]
head(std_ibd_df_final)

#Let's make it linear model ready
dim(std_ibd_df_final) # 97962     3

# Long format to wide
std_ibd_WIDE <- spread(std_ibd_df_final, MARKER, std_ibd)

# Final table
dim(std_ibd_WIDE)   # 519 190
head(std_ibd_WIDE)

# Save std_ibd_WIDE
saveRDS(std_ibd_WIDE, file = "./results/std_IBD_WIDE.Rds")

# Session Info
writeLines(capture.output(sessionInfo()), "./results/sessionInfo.txt")
