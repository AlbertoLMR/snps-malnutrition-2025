# Title "Binding two databases with IDs and sex IDs (without duplicates)"

# Author: Alberto López Martínez Rojas
# Contact: albertolomr@06gmail.com

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

# After running this script I compared the new CSV with the old .fam to validate how many
# Samples I had. 
