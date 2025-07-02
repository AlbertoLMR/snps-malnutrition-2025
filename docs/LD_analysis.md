# Linkage Disequilibrium Analysis Documentation

## Overview

This document provides detailed documentation for conducting a Linkage Disequilibrium (LD) analysis using R and PLINK. The analysis includes data pruning, complete dataset analysis, sex-stratified analysis, and a SNP summary table generation. 

## Dependencies 

The dependencies used for the analysis and correction are listed below: 

### Software

- R programming language for Windows: [Download R-4.5.1 for Windows.](https://cran.r-project.org/bin/windows/base/)
- RStudio Desktop: [RStudio Desktop - Posit](https://posit.co/download/rstudio-desktop/)
- Excel or LibreOffice Calc: [Microsoft 365](https://www.office.com/), [Download LibreOffice](https://www.libreoffice.org/download/download-libreoffice/)

### R packages

- R dependency: `tidyverse`
- PLINK: [PLINK 1.9](https://www.cog-genomics.org/plink/)

## Input files required 

The LD analysis requires three PLINK binary files with the same prefix:

- `.bed` file: Binary genotype data

- `.bim` file: SNP information (chromosome, position, alleles)

- `.fam` file: Sample information including sex coding (1=male, 2=female)

## Workflow overview

The LD analysis provided in `code/LD_analysis_example.R` follows a structured workflow: 

1. **Data preparation**: Set up workspace and create output directories
2. **SNP pruning**: Remove highly correlated SNPs using pairwise LD
3. **Complete Dataset Analysis**: Analyze LD across all samples
4. **Sex-Stratified Analysis**: Separate analysis for males and females
5. **Data Filtering**: Filter results by LD strength categories
6. **SNP Summary**: Generate allele information table

## Detailed step-by-step guide

### Step 1: Workspace Preparation

```R
# Load required packages
library(tidyverse)

# Set working directory to your project folder
setwd("path/to/your/directory")

# Create results directory
if (!dir.exists("LD_results")) {
  dir.create("LD_results")
}
```

**Purpose:** Ensures your R environment is properly configured and output files are organized. 

### Step 2: SNP pruning

```R
system(paste0("./plink.exe --bfile input_prefix", 
              " --indep-pairwise 50 5 0.5",
              " --out output_prefix"))
```

**Parameters Explained**:

`--indep-pairwise 50 5 0.5`: Pruning algorithm settings

- `50`: Window size (number of SNPs to consider at once)
- `5`: Step size (shift window by this many SNPs each iteration)
- `0.5`: R² threshold (remove SNPs with R² > 0.5)

**Output**: Creates a `.prune.in` file containing SNPs that pass the pruning filter.

**Why Pruning?**: Removes highly correlated SNPs that could bias downstream analyses while retaining representative SNPs across the genome.

### Step 3: Create pruned dataset

```R
system(paste0("./plink.exe --bfile input_prefix",
              " --extract file_name.prune.in",
              " --make-bed --out pruned_prefix"))
```

**Purpose**: Creates new binary files containing only the pruned SNPs for LD analysis.

### Step 4: LD analysis - complete dataset

```R
system(paste0("./plink.exe --bfile pruned_prefix", 
              " --r2 --ld-window-kb 1000 --ld-window r2 0.05",
              " --out LD_results/complete_ld"))
```

**Parameters Explained**:

- `--r2`: Calculate r² values for LD measurement
- `--ld-window-kb 1000`: Consider SNP pairs within 1000 kb distance
- `--ld-window-r2 0.05`: Only report SNP pairs with r² ≥ 0.05

**Output**: `.ld` file containing LD results. Example LD table: 

| CHR_A | BP_A | SNP_A | CHR_B | BP_B | SNP_B |  R2  |
| :---: | :--: | :---: | :---: | :--: | :---: | :--: |
|   1   | ###  |  rs#  |   1   | ###  |  rs#  | 0.06 |

- **CHR_A**: Chromosome of first SNP
- **BP_A**: Base position of first SNP
- **SNP_A**: Identifier of first SNP
- **CHR_B**: Chromosome of second SNP
- **BP_B**: Base position of second SNP
- **SNP_B**: Identifier of second SNP
- **R2**: Linkage disequilibrium r² value (0-1)

### Step 5: Data processing and filtering

The script processes LD results and filters them into meaningful categories:

**High LD**: r² between 0.3 and 0.5

- Indicates moderate to strong linkage

**Low LD**: r² between 0.05 and 0.3

- Indicates weak to moderate linkage

```R
# Example filtering
high_ld <- ld_data %>% filter(R2 >= 0.3 & R2 <= 0.5)
low_ld <- ld_data %>% filter(R2 >= 0.05 & R2 < 0.3)
```

### Step 6: Sex-stratified analysis

```R
# Males only
system(paste0("./plink.exe --bfile pruned_prefix",
              " --filter-males --r2 --ld-window-kb 1000 --ld-window-r2 0.05",
              " --out LD_results/males_ld"))
```

```R
# Females only
system(paste0("./plink.exe --bfile pruned_prefix",
              " --filter-females --r2 --ld-window-kb 1000 --ld-window-r2 0.05",
              " --out LD_results/females_ld"))
```

**Purpose**: Analyzes LD patterns separately by sex, which can reveal sex-specific genetic architecture differences.

### SNP Summary table

The final section creates a comprehensive summary of SNPs that passed pruning:

```R
# Read pruned .bim file to get allele information for the SNPS
# Note: Use the same prefix as your pruned dataset from step 2
bim_data <- read.table("file_name_input_prefix.bim", header = FALSE)
colnames(bim_data) <-c("CHR", "SNP", "CM", "BP", "A1", "A2") #add standard .bim column names

# Explanation of .bim file columns:
# CHR = Chromosome number
# SNP = SNP identifier 
# CM = Genetic distance (centimorgans)
# BP = Base-pair position
# A1 = Allele 1 (usually minor allele)
# A2 = Allele 2 (usually major allele)

#create SNP summary table 
snps_summary <- bim_data %>%
  mutate(
    Allelic_Variants = 2, 
    Base_Pairs = paste(A1, A2, sep = "/") #format: A/T, G/C, etc
  ) %>%
  
  select(SNP_Name = SNP, Allelic_Variants, Base_Pairs) %>%
  arrange(SNP_Name)
```

**Use Cases**:

- Quality control verification
- Allele frequency analysis preparation
- Documentation 

## Customization options

### Adjusting pruning parameters

- **More stringent pruning**: Decrease R² threshold (e.g., 0.2)
- **Less stringent pruning**: Increase R² threshold (e.g., 0.8)
- **Larger windows**: Increase window size for broader LD assessment

### Modifying LD analysis parameters

- **Distance limits**: Adjust `--ld-window-kb` for different distance ranges
- **Reporting threshold**: Modify `--ld-window-r2` to change minimum reporting threshold

### Custom filtering categories

Users can define their own r² ranges based on study requirements:

```
very_high_ld <- ld_data %>% filter(R2 >= 0.8)
custom_range <- ld_data %>% filter(R2 >= 0.1 & R2 <= 0.4)
```

## Troubleshooting

### Common isues

**PLINK not found**: Ensure PLINK executable is in your working directory or system PATH.

**Empty results**: Check that:

- Input files are properly formatted
- Sex coding is correct in .fam files (1 = male, 2 = female)
- Pruning didn't remove all SNPs

## Output organization

The script generates organized outputs in the `LD_results/` directory:

- Raw LD results (.ld files)
- Processed CSV files for easy viewing
- Filtered results by LD strength
- SNP summary table

## Best practices

**1. Backup original data** before running analysis

**2. Document parameter choices** and rationale

**3. Validate results** using known LD patterns if available

## Applications

This LD analysis is commonly used for:

- **GWAS preparation**: Identifying independent SNPs
- **Population genetics**: Understanding genetic structure