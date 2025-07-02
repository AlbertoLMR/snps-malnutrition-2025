# Fam file Sex ID Update Documentation

## Overview

This document provides detailed documentation for updating PLINK .fam files with incorrect sex ID coding. The workflow addresses the common issue where .fam files contain 0 values for sex identification instead of the required PLINK format (1 = male, 2 = female). The analysis includes database merging, ID matching validation, sex code updating, and creation of clean binary datasets ready for genetic analysis.

## Dependencies

The dependencies used for the analysis and correction are listed below:

### Software

- R programming language for Windows: [Download R-4.5.1 for Windows.](https://cran.r-project.org/bin/windows/base/)
- RStudio Desktop: [RStudio Desktop - Posit](https://posit.co/download/rstudio-desktop/)
- Excel or LibreOffice Calc: [Microsoft 365](https://www.office.com/), [Download LibreOffice](https://www.libreoffice.org/download/download-libreoffice/)

### R packages

No additional R packages required beyond base R functionality.

## Input files required

The FAM file update process requires the following input files:

### Primary Input Files

- Original PLINK binary files

   (same prefix):

  - `.bed` file: Binary genotype data
  - `.bim` file: SNP information
  - `.fam` file: Sample information with incorrect sex coding (contains 0 values)

### Reference database files

In my case I had 2 different databases with the SNPs IDs.

- **sex_sheet1.csv**: First database containing individual IDs and correct sex information

- **sex_sheet2.csv**: Second database containing individual IDs and correct sex information

### Expected file format

Reference CSV files should contain at minimum:

- **ID column**: Individual identifiers matching those in [.fam file](https://genomicsbootcamp.github.io/book/genotype-files-in-practice.html)	
- **Sex column**: Correct sex coding (1=male, 2=female)

## Workflow overview

The FAM file sex ID update workflow provided in the scripts follows a structured approach:

1. **Database Preparation**: Merge multiple sex reference databases and remove duplicates
2. **ID Validation**: Check matching rates between reference data and original .fam file
3. **Sex Code Integration**: Use spreadsheet software to match and update sex information
4. **Data Filtering**: Remove individuals with unknown sex coding
5. **Binary File Creation**: Generate updated PLINK binary files with correct sex information

## Detailed step-by-step guide

### Step 1: Database preparation

```R
#clear work space
rm(list = ls()) 

# Check current directory and set working directory
getwd()
setwd("path/to/your/directory")       

# Read first CSV file
sheet1 <- read.csv("sex_sheet1.csv")
# Read second CSV file
sheet2 <- read.csv("sex_sheet2.csv")

# Combine both datasets by rows
combined <- rbind(sheet1, sheet2)

# Remove duplicates based on first column (ID), keeping first occurrence
combined_unique <- combined[!duplicated(combined[,1]), ]

# Print row counts for each step
print(paste("Sheet 1 had:", nrow(sheet1), "rows"))
print(paste("Sheet 2 had:", nrow(sheet2), "rows"))
print(paste("Combined:", nrow(combined), "rows"))
print(paste("After removing duplicates:", nrow(combined_unique), "rows"))

# Save cleaned dataset as CSV
write.csv(combined_unique, "sex_reference_combined.csv", row.names = FALSE)
```

**Purpose:** Creates a comprehensive, duplicate-free reference database containing correct sex information from multiple sources.

**Quality Check:** The script reports:

- Original row counts from each database
- Combined total before deduplication
- Final count after duplicate removal
- Number of duplicates eliminated

### Step 2: ID Matching validation

```R
#clear work space
rm(list = ls()) 

# Check current directory and set working directory
getwd()
setwd("path/to/your/directory")  

# Read updated CSV and original .fam 
sex_reference <- read.csv("sex_reference_combined.csv")
fam_data <- read.csv("old_fam.csv") 

# Extract ID columns from both datasets
sex_ids <- sex_reference$"iid"
fam_ids <- fam_data$"iid"

# Find IDs that exist in both datasets
common_ids <- intersect(sex_ids, fam_ids)
# Print match statistics
print(paste("Exact matches found:", length(common_ids)))
print(paste("Sex reference file has", length(sex_ids), "IDs"))
print(paste("FAM data file has", length(fam_ids), "IDs"))

# Calculate percentage of FAM IDs that have matches
print(paste("Match percentage:", round(length(common_ids)/length(fam_ids)*100, 1), "%"))

# Evaluate match quality and provide recommendations
print("=== SUMMARY ===")
if(length(common_ids)/length(fam_ids) >= 0.8) {
  print("Most IDs match. Ready to proceed with merging.")
} else if(length(common_ids)/length(fam_ids) >= 0.7) {
  print("GOOD match rate, but some IDs don't match. Check the mismatches above.")
} else {
  print("low matchtch rate. Need to fix ID formatting before proceeding.")
}

print(paste("Ready for merging:", length(common_ids), "samples"))
```

**Match Quality Assessment:**

- **≥80% match rate**: Excellent - ready to proceed
- **70-79% match rate**: Good - some manual ID checking recommended
- **<70% match rate**: Poor - ID formatting issues need resolution

**Purpose:** Validates that sufficient overlap exists between reference data and .fam file before proceeding with sex code updates.

### Step 3: Sex code integration

This step uses spreadsheet software (LibreOffice Calc or Excel) to perform the actual sex code matching between the original .fam and the merged database:

**LibreOffice Calc Formula:**

`=IFERROR(VLOOKUP(B2;$I$2:$J$800;2;0);E2 #using the columns on my spreadsheet`

**Formula Explanation:**

- `VLOOKUP(B2;$I$2:$J$800;2;0)`: Looks up ID in B2 within reference range I2:J800, returns sex value from column 2
- `IFERROR(...;E2)`: If no match found, keeps original value from column E2
- `$I$2:$J$800`: Absolute reference to sex database range (adjust range as needed)

**Process:**

1. **Import data**: Load both .fam file and sex reference database into spreadsheet
2. **Apply formula**: Use VLOOKUP to match IDs and update sex codes
3. **Preserve order**: Maintains original sample order and Family IDs
4. **Data integrity**: Only updates records with exact ID matches
5. **Export results**: Save as "updated_corrected_FAM.csv"

### Step 4: Data filtering and cleanup

```R
#clear workspace
rm(list = ls()) 

# Check current directory and set working directory
getwd()
setwd("path/to/your/directory")  

# Read the manually updated .fam file with corrected sex IDs
fam_data <- read.csv("updated_corrected_FAM.csv")
# Remove rows where sex_id is 0 (unknown/uncoded)
fam_filtered <- fam_data[fam_data$sex_id != 0, ]

# Print filtering results
print(paste("Filtered rows:", nrow(fam_filtered)))
print("Sex distribution:")
# Show count of males (1) and females (2)
print(table(fam_filtered$sex_id))

# Save as .fam format (tab-separated, no headers)
write.table(fam_filtered, "final.fam", quote = FALSE,  row.names = FALSE, col.names = FALSE, sep = "\t")
```

**Purpose:** Removes individuals with incomplete sex information and creates a clean .`fam` file in proper PLINK format.

**Quality Control:** Reports final sex distribution to verify expected male/female ratios.

### Step 5: PLINK Binary File Creation

```R
#clear workspace
rm(list = ls()) 

# Check current directory and set working directory
getwd()
setwd("path/to/your/directory")  

# Read the final .fam file without headers
fam_corrected <- read.table("final.fam", header = FALSE)
# Create sex update file with Family ID, Individual ID, and Sex columns
sex_update <- fam_corrected[, c(1, 2, 5)]  # FID, IID, SEX columns of a standard .fam
# Save sex update file for PLINK
write.table(sex_update, "sex_update.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Create keep list with individuals who have sex coded (all remaining individuals)
keep_list <- fam_corrected[fam_corrected$V5 != 0, c(1,2)]

# Save keep list for PLINK filtering
write.table(keep_list, "keep_sex_coded.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Set input and output file prefixes for PLINK
input_prefix <- "old"
output_prefix <- "final"

# Run PLINK to update sex information and create new binary files
system(paste0("./plink.exe --bfile ", input_prefix,
              " --update-sex sex_update.txt",
              " --keep keep_sex_coded.txt",
              " --make-bed --out ", output_prefix))
```

**PLINK Parameters Explained:**

- `--bfile old`: Input binary files with prefix "old"
- `--update-sex sex_update.txt`: Updates sex information from specified file
- `--keep keep_sex_coded.txt`: Retains only individuals listed in keep file
- `--make-bed`: Creates new binary files (.bed/.bim/.fam)
- `--out final`: Output files with prefix "final"

**Output:** Creates updated PLINK binary files ready for genetic analysis.

## File Format Specifications

### Standard .fam File Format

PLINK `.fam` files require 6 columns (tab-separated, no header):

| Column |   Name    |  Description  | Example |
| :----: | :-------: | :-----------: | :-----: |
|   1    |    FID    |   Family ID   | FAM001  |
|   2    |    IID    | Individual ID | IND001  |
|   3    |    PID    |  Paternal ID  |    0    |
|   4    |    MID    |  Maternal ID  |    0    |
|   5    |    SEX    |   Sex code    | 1 or 2  |
|   6    | PHENOTYPE |   Phenotype   |   -9    |

### Sex Coding Standards

- **1**: Male
- **2**: Female
- **0**: Unknown (problematic for most analyses)
- **-9**: Missing (alternative unknown coding)

## Troubleshooting

### Common Issues

**Low ID Match Rates**

- Check ID formatting consistency between files
- Verify column names match exactly in R code
- Look for leading/trailing spaces in IDs
- Check for case sensitivity differences

**VLOOKUP Not Working in Spreadsheet**

- Ensure reference range is correctly specified with absolute references ($)
- Check that ID columns are formatted consistently (text vs. numeric)
- Verify lookup range contains the expected data

**PLINK Errors**

- Ensure PLINK executable is in working directory or system PATH
- Check that input binary files exist and are readable
- Verify file paths are correct in system commands

**Empty Results After Filtering**

- Check that sex_id column name matches exactly in filtering command
- Verify sex codes are numeric (1, 2) not text ("1", "2")
- Ensure updated `.fam` file was saved correctly

## Output Organization

The workflow generates several key output files:

### Intermediate Files

- `sex_reference_combined.csv`: Merged reference database
- `updated_corrected_FAM.csv`: `.fam` file with updated sex codes
- `sex_update.txt`: PLINK sex update file
- `keep_sex_coded.txt`: PLINK keep list

### Final Output Files

- `final.fam`: Clean file with proper sex coding
- `final.bed`: Updated binary genotype file
- `final.bim`: SNP information file (unchanged)

## Best Practices

**1. Backup original files** before starting any modifications

**2. Document ID matching process** and any manual corrections made

**3. Use sex distributions** that make biological sense for your population

**4. Test with small subset** before processing large datasets

**5. Keep detailed logs** of match rates and filtering decisions

## Applications

This `.fam` file update process is essential for:

- **GWAS preparation**: Ensuring proper sex coding for analysis
- **Linkage analysis**: Family-based genetic studies
- **Quality control**: Data cleaning and standardization

## Quality Control Recommendations

### Pre-Processing Checks

- Verify reference databases contain expected sex ratios
- Check for systematic ID formatting differences
- Validate that all required individuals have sex information

### Post-Processing Validation

- Compare final sex distribution to expected population ratios
- Verify sample sizes match expectations after filtering
- Test updated files with downstream analysis tools

### Documentation Requirements

- Record original and final sample counts
- Document match rates and filtering decisions
- Note any manual corrections or exceptions made
