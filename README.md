# DepMap_correlation_tools (public version)
by Heyuan Li

This is a correlation tools to analyse and visualize DepMap data, in order to discover novel genetic interactions and potential therapeutic targets/biomarkers. While most of the functions are also available on the DepMap.org portal, this tool highlight the functions to search for tissue specific correlations, and correlations for a list of genes. 

## instructions

After downloading the data from DepMap.org, run `data_formatting.R` to to generate global environment. 
To use the functions, run `sourse('functions.R')`. 

## Data used

All the data are organized as cell line in rows, genes/pathways in columns. Here I used CRISPR screen dependency scores for most of the analysis for illustration. It can be swapped to other data tables:

### Cell line & genes
`ccle_cn`: copy number variations in CCLE
`ccle_rnaseq`: RNA expression in CCLE
`ccle_protein`: Protein expression in CCLE
`drive_rsa`: RSA score in project DRIVE (Robert et al, Cell, 2017)
`shrna_d2`: DEMETER2 score in combined shRNA screens

### Cell line & pathways
`gsva`: gsva score calculated from `ccle_rnaseq`


### data sources
`drive_rsa` is a replicate analysis from the original paper (Robert et al, Cell, 2017)
`gsva` is calculated from `gsva()` from GSVA package from Bioconductor

The rest of these data are downloaded from DepMap.org (22Q1) and are pre-processed by `data_formmating.R`
