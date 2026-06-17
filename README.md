# Genetic modifiers of microbiome traits in IBD families

Analysis code for the study *Linkage analysis identifies novel genetic modifiers of microbiome traits in families with inflammatory bowel disease* (Sharma et al., *Gut Microbes*, 2022).

Paper: https://doi.org/10.1080/19490976.2021.2024415 (open access, CC BY-NC 4.0)

## Overview

Using SNP genotypes together with 16S rRNA gene microbiome and phenotype data from the Kiel IBD kindred cohort (IBD-KC), this study examined the relationship between genetic and microbiome similarity in families of IBD patients, and then performed a genome-wide quantitative trait locus (QTL) linkage analysis of microbiome traits. The code here covers the four analysis stages: pedigree-based heritability of genus abundances, pairwise genetic and microbiome similarity at IBD risk loci, the genome-wide QTL linkage scan, and SNP-wise association within the linked regions.

## Repository structure

The analysis scripts are in `R/`, numbered in the order they are run:

- `01_heritability.R`: pedigree-based heritability of genus abundances (MERLIN variance components)
- `02_ibd_calculation.R`: identity-by-descent sharing at the IBD risk loci (MERLIN)
- `03_pairwise_similarity.R`: pairwise genetic and microbiome similarity among family pairs
- `04_linkage_analysis.R`: genome-wide QTL linkage scan (MERLIN-REGRESS)
- `05_extract_interval.R`: two-unit support intervals around the linkage peaks
- `06_Association_2USI.R`: SNP-wise association within the linkage regions (WISARD/GEMMA)
- `07_locus_zoom_plots.R`: locus-zoom plots of LOD score, association, and genes
- `08_linkage_mGWAS_overlap.R`: overlap of the linkage regions with previously reported microbiome-GWAS SNPs

## Data

The individual-level genotype and microbiome data are not included, as they cannot be shared for data protection reasons. The scripts are provided to document the analysis rather than to reproduce it end to end, so they will not run as-is without the underlying IBD-KC data.

## Tools and packages

The analysis was carried out in R. Packages used across the scripts include dplyr, tidyr, ggplot2, patchwork, forcats, data.table, pedtools, and biomaRt (v2.44.4). The genome-wide linkage scan was run in MERLIN (merlin-regress, pedstats), the within-region SNP association in WISARD with GEMMA, and SNP positions for the mGWAS comparison were obtained from Ensembl VEP. Add the R and MERLIN version numbers and the package versions, for example from the output of sessionInfo()

## Citation

Sharma A, Szymczak S, Rühlemann M, Freitag-Wolf S, Knecht C, Enderle J, Schreiber S, Franke A, Lieb W, Krawczak M, Dempfle A. Linkage analysis identifies novel genetic modifiers of microbiome traits in families with inflammatory bowel disease. *Gut Microbes*. 2022;14(1):2024415. doi:10.1080/19490976.2021.2024415. PMID: 35129060.

```bibtex
@article{sharma2022linkage,
  title   = {Linkage analysis identifies novel genetic modifiers of microbiome traits in families with inflammatory bowel disease},
  author  = {Sharma, Arunabh and Szymczak, Silke and R{\"u}hlemann, Malte and Freitag-Wolf, Sandra and Knecht, Carolin and Enderle, Janna and Schreiber, Stefan and Franke, Andre and Lieb, Wolfgang and Krawczak, Michael and Dempfle, Astrid},
  journal = {Gut Microbes},
  volume  = {14},
  number  = {1},
  pages   = {2024415},
  year    = {2022},
  doi     = {10.1080/19490976.2021.2024415},
  pmid    = {35129060}
}
```


