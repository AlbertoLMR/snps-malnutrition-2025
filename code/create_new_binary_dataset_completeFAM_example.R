# Title "Create new binary dataset for the updated and correct .fam"

# Author: Alberto López Martínez Rojas
# Contact: albertolomr@06gmail.com

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