# Title "Delete sex ID values with 0 in new updated .fam"

# Author: Alberto López Martínez Rojas
# Contact: albertolomr@06gmail.com

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

# The next step is to create a new binary set for the updated and complete .fam
