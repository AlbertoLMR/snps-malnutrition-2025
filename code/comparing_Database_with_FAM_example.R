# Title "Check if IDs match between databases and old .fam data "

# Author: Alberto López Martínez Rojas
# Contact: albertolomr@06gmail.com

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

# After running this script I used Excel/Libre Calc to match the correct Sex Id from the database
# within the old .fam. Creating a new updated .fam (still with some SNPs with a value of 0 for 
# Sex ID)

# After this step I added the standard column names for a .fam file manually on Excel/Libre Calc
# to this new CSV and the old.fam CSV