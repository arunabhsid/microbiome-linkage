################################################################################

#     LOCUS-ZOOM PLOT FOR LINKAGE & FAMILY BASED ASSOCIATION ANALYSIS

################################################################################

# Locus-zoom plotting approach adapted from a tutorial by Sina Rüeger
# Requires biomaRt (Bioconductor):
   # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   # BiocManager::install("biomaRt")

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
}

# libraries
library(dplyr) 
library(forcats) #factors
library(ggplot2)
library(patchwork) 
library(biomaRt)

## Region to plot (edit these three per region)
region_file  <- "./results/two_unit_interval/odoribacterCHR22_2USI_region1.txt"
assoc_file   <- "./results/two_unit_interval/assoc_pvals/odoribacter_gemma_2USI_region1.txt"
region_label <- "Odoribacter_CHR22_region1"
highlight_bp <- numeric(0)   # IBD-associated SNP(s) to mark red for this region; numeric(0) if none
# regions that have one:
#   barnesiella chr4          : c(102865304, 103434253)
#   clostridium chr14 region1 : 88472595
#   roseburia chr19           : 55380214
#   odoribacter chr22 region2 : 33229830

################################################################################

#                     Two-unit support interval region plot

################################################################################

# load the data
# trait_name is the name of df with two-unit support interval (2USI)

# for other traits: sep = "\t"
trait_name <- read.table(region_file, 
                         sep = "\t",
                         header = TRUE)

# specifying start and end positions in bp
start_position_in_BP <- trait_name$BP[which.min(trait_name$BP)]
end_position_in_BP <- trait_name$BP[which.max(trait_name$BP)]

# table
start_position_in_cM <- trait_name$cM[which.min(trait_name$cM)]
end_position_in_cM <- trait_name$cM[which.max(trait_name$cM)]
size_of_region_cM <- (end_position_in_cM - start_position_in_cM)

# if the chromosomal position is in basepairs, convert it to megabases
trait_name$MB <- trait_name$BP / 1000000
start_position_in_MB <- trait_name$MB[which.min(trait_name$MB)]
end_position_in_MB <- trait_name$MB[which.max(trait_name$MB)]

size_of_region <- (end_position_in_MB - start_position_in_MB)

# extract information regarding the chromosome for the plot and biomaRt
chr_number <- trait_name$CHR[which.max(trait_name$CHR)]

# get a name of the phenotype for the plot
trait_name$PHENOTYPE <- gsub(" ", "", trait_name$PHENOTYPE)
phenotype_name <- trait_name$PHENOTYPE[1]

# p1 is a 2USI plot.
# intercept for the significance level
p1 <- ggplot(data = trait_name) + geom_line(aes(MB,LOD)) + theme_minimal() +
  labs(subtitle = paste("Two-unit support interval for", phenotype_name, "on chromosome", format((chr_number)), "from", format((start_position_in_MB)), "to", format((end_position_in_MB)), "Mb")) 

p1 <- p1 + geom_hline(yintercept = 3, linetype = "dashed", color = "red")

# arrow marking the interval, dashed drop-line at the peak (derived from the loaded region)
arrow_y   <- 1
peak_mb   <- trait_name$MB[which.max(trait_name$LOD)]
peak_lod  <- max(trait_name$LOD)
p1 <- p1 +
  geom_segment(aes(x = start_position_in_MB, xend = end_position_in_MB, y = arrow_y, yend = arrow_y),
               arrow = arrow(length = unit(0.3, "cm"), ends = "both"), size = 0.5) +
  geom_segment(aes(x = peak_mb, xend = peak_mb, y = peak_lod, yend = arrow_y), linetype = 2, size = 0.2)

# Print
print(p1)

################################################################################

#                          Gene annotation plot

################################################################################

#extract gene annotations using biomaRt
#select biomaRt database and dataset hosted by Ensembl first

gene_ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)

################################################################################
#                     gene annotations in the region (all biotypes)
################################################################################

output_genes_region <- getBM(
  attributes = c('start_position', 'end_position', 'ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
  filters    = c('chromosome_name', 'start', 'end'),
  values     = list(chr_number, start_position_in_BP, end_position_in_BP),
  mart       = gene_ensembl)

# remove duplicate gene names
output_genes_region <- output_genes_region[!duplicated(output_genes_region$external_gene_name), ]
################################################################################

# REORDER genes
output_genes_region <- output_genes_region %>% mutate(external_gene_name = fct_reorder(external_gene_name,
                                                                                       start_position, .desc = TRUE))
# for paper: write to file
write.csv(output_genes_region, file = paste0("./results/gene_list/", region_label, "_genes.csv"), quote = F, row.names = F)

# add position in Mb column
output_genes_region$MBstart_position <- output_genes_region$start_position / 1000000
output_genes_region$MBend_position <- output_genes_region$end_position / 1000000

# create a gene annotation plot
ggplot(data = output_genes_region) +
  geom_linerange(aes(x = external_gene_name, ymin = MBstart_position, ymax = MBend_position)) +
  coord_flip() +
  ylab("")

# define the plot range for the x-axis
plot_range <- c(min(start_position_in_MB, output_genes_region$MBstart_position),
                max(end_position_in_MB, output_genes_region$MBend_position))

# group protein coding genes and color them
output_genes_region <- output_genes_region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype),
                                                                                     "protein_coding"), external_gene_name = fct_reorder2(external_gene_name,
                                                                                                                                          MBstart_position, gene_biotype_fac, .desc = TRUE))
# rename x-axis as "genes"
output_genes_region1 <- output_genes_region %>% rename(
  genes = external_gene_name
)

# gene annotation plot (without gene names - just lines)
p2 <- ggplot(data = output_genes_region1) +
  geom_linerange(aes(x = genes, ymin = MBstart_position, ymax = MBend_position)) +
  coord_flip() + xlab("Genes") + ylab("Mb") + 
  ylim(plot_range) +
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "",
        panel.grid.major.y = element_blank()) +
  expand_limits(y=c(-1,1))

print(p2)
################################################################################

#                          Association plot

################################################################################

# assoc_data is a summary statistics study 
assoc_data <- read.table(assoc_file, 
                         sep = "\t",
                         header = TRUE)
#rename columns 
assoc_data <- assoc_data %>% 
  rename(
    BP = POS,
    P = P_LRT)

# specifying start and end positions in bp
start_position_in_BP <- assoc_data$BP[which.min(assoc_data$BP)]
end_position_in_BP <- assoc_data$BP[which.max(assoc_data$BP)]

# if the chromosomal position is in basepairs only in your data, you can convert it into megabases
assoc_data$MB <- assoc_data$BP / 1000000
start_position_in_MB <- assoc_data$MB[which.min(assoc_data$MB)]
end_position_in_MB <- assoc_data$MB[which.max(assoc_data$MB)]

# p3 is a plot for summary statistics 
p3 <- ggplot(data = assoc_data) + geom_point(aes(MB,-log10(P))) +theme_minimal()  
print(p3)

################################################################################
#                   highlight SNPs of interest
################################################################################

# SNPs to mark are set via highlight_bp in the header (numeric(0) = none)
highlight_df <- assoc_data %>% filter(BP %in% highlight_bp)

# summary statistics, with SNPs of interest marked in red
p3 <- ggplot(data = assoc_data) + geom_point(aes(MB, -log10(P))) +
  geom_point(data = highlight_df, aes(MB, -log10(P)), color = "red", size = 1) +
  theme_minimal()
print(p3)

################################################################################

#                          Final plot

################################################################################

# combine three plots
jpeg(paste0("./plots/locus_zoom/", region_label, ".jpg"),height=4,width=10,units='in',res=600)
p1b <- p1 + xlab("")+ theme(axis.title.x=element_blank(), axis.text.x = element_blank()) 
p3b <- p3 + xlab("")+ theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + xlim(plot_range) 
p1b + p3b + p2 + plot_layout(nrow = 3, heights = c(1, 1, 0.5))
dev.off()