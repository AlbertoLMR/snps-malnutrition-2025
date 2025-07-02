# Linkage Disequilibrium analysis for SNPs

This repository provides tools for conducting Linkage Disequilibrium (LD) analysis on human SNP arrays using R and PLINK, with utilities for correcting common `.fam` file issues

## Overview

The goal of the project is to analyze malnutrition related SNPs through LD analysis using R scripts and PLINK specific R dependencies, as well as implementing a sex-stratified LD analysis prior to Genome-wide association studies (GWAS). 

A secondary objective is the use of R scripts and other tools to solve a problem related to the correction of `.fam files` (file extension necessary to the LD analysis via PLINK). 

## Features

- LD analysis for the whole dataset with R and PLINK
- Sex-stratified LD analysis with R and Plink
- `.fam file` sex coding correction and dataset cleanup

## Dependencies

The dependencies used for the analysis and correction are listed below: 

### Software

- R programming language for Windows: [Download R-4.5.1 for Windows.](https://cran.r-project.org/bin/windows/base/)
- RStudio Desktop: [RStudio Desktop - Posit](https://posit.co/download/rstudio-desktop/)
- Excel or LibreOffice Calc: [Microsoft 365](https://www.office.com/), [Download LibreOffice](https://www.libreoffice.org/download/download-libreoffice/)

### R packages

- R dependency: `tidyverse`
- PLINK: [PLINK 1.9](https://www.cog-genomics.org/plink/)

## LD analysis using PLINK

To use PLINK for LD analysis you need to have three extension files: 

1. `.fam` file
2. `.bim` file
3. `.bed`file 

PLINK uses this basic command to compute the LD analysis and output a result table

```R
system(paste0("PLINK --bfile file_name_input_prefix", " --r2 --ld-window-kb --ld-window-r2"," --out file_name_output_prefix"))
```

In my case, I installed PLINK on my project directory and I was able to run it like: 

```R
system(paste0("./plink.exe --bfile file_name_input_prefix", " --r2 --ld-window-kb --ld-window-r2"," --out file_name_output_prefix"))
```

## .fam file sex coding issue 

### Problem

PLINK `.fam` files can contain incorrect sex coding (0 for unknown/missing) instead of the standard format (1=male, 2=female). This can cause issues with downstream analyses that require proper sex information. 

### Solution approach

This repository includes utilities to: 

- Cross-reference `.fam` files with external databases containing correct sex information
- Update sex coding from 0 to 1/2 values based on SNPs ID matching
- Remove samples not found in the reference database or with missing sex coding
- Generate corresponding `.bed`and `.bim` files to maintain file consistency

### Case example

My specific implementation addresses a scenario where:

- Original `.fam` file had all sex values coded as 0 
- External SNP database contained samples with correct sex coding (1/2)
- Solution: match samples by ID, update sex coding, and create clean dataset for LD analysis

See `docs/FAM_file_sex_ID_update.md` for detailed implementation of the solution scripts and utilities. 

## Repository structure

- The directory `code/` contains the analysis scripts as well as the solution scripts
- The directory `docs/` contains the detailed documentation for the implementation of the scripts and utilities 

## Contact

- email: albertolomr06@gmail.com

## License

