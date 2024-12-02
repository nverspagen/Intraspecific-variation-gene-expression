# Intraspecific-variation-gene-expression
Code and metadata used for data analysis for the manuscript: "Mild warming induces divergent plastic responses in gene expression among populations of a temperate butterfly"

# Code
Code consists of four files:
- DataExploration.Rmd. This script performs the data exploration and quality check. It needs the raw sequencing data that are cleaned and aligned to the reference genome (https://doi.org/10.5524/100915) through https://nf-co.re/rnaseq/3.11.2 to run.
- WGCNAAnalysis.R. This script does some more data exploration, and performs WGCNA analysis. It needs output from DataExploration.Rmd to run.
- DifferentialExpression.Rmd. This script performs differential expression analysis. It needs output from DataExploration.Rmd to run.
- This script performs enrichment analysis. It needs outputs from WGCNAAnalysis.R and DifferentialExpression.Rmd, as well as reference genome annotation files (https://doi.org/10.5524/100915) to run.

# Data
Small datafiles are included here.
- Data_matrix_RNAseq.csv includes all information on the samples. ID column provides a unique identifier for each individual. It is made up of population, family, temperature, clutch, replicate and individual. The population codes for country of origin, with AF for Finland and CAT for Spain. The Family column indicates the ID number of the mother, which is our definition of an experimental family. Temperature shows day rearing temperature. Clutch shows the egg clutch ID and replicate further defines the experimental group the individual grew up in. Sex shows the sex of the individual as determined through SNPs on the Z chromosome. Other columns provide further information on the individual, but were not used in further analysis.
- F22FTSEUHT0133_MELdeyuR.xlsx provides information on sequencing such as library ID, flow cell name and lane number. Sample name corresponds to the sample ID from the data matrix.
Data matrix and sequencing info are combined and used for data exploration. Temperature and population identifiers are also used in further analysis.
- blast_results_final.txt are the results of blast against drosophila genome.
