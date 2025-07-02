# Title Linkage Disequilibrium analysis with PLINK

# Author: Alberto López Martínez Rojas
# Contact: albertolomr@06gmail.com

# Dependencies
  # tidyverse
  # PLINK

#-----------------------Prepare work space for LD analysis----------------------

#install and load packages
install.packages("tidyverse") #one time only
library(tidyverse) 

#clear workspace
rm(list = ls()) 

#get the current working directory 
getwd()

#set the correct working directory 
setwd("path/to/your/directory")       

#create LD results folder to save results (optional)
if (!dir.exists("LD_results")) {
    dir.create("LD_results")
  } #the if statement is to always ensure a results directory exists

#------------------------LD analysis for complete data set----------------------

#first step: pruning data to remove highly correlated SNPs (through pairwise)
system(paste0("./plink.exe --bfile file_name_input_prefix", 
              " --indep-pairwise 50 5 0.5",
              " --out file_name_output_prefix"))

# 50 =  window size (50 SNPs)
# 5 = shift by 5 SNPs each iteration of the analysis
# 0.5 = keep SNPs that have a value of R2 <= 0.5

#this PLINK command outputs a file with a .prune.in file extension

#second step: create the pruned data files - .bed, .bim and, .fam with the pruned SNPs
system(paste0("./plink.exe --bfile file_name_input_prefix",
              " --extract file_name.prune.in",
              " --make-bed --out file_name_output_prefix"))

#this PLINK command outputs a new pruned binary data set to use for the LD analysis

#third step: LD analysis for the complete binary pruned data set
system(paste0("./plink.exe --bfile file_name_input_prefix",
              " --r2 --ld-window-kb 1000 --ld-window-r2 0.05",
              " --out LD_results/file_name_output_prefix"))

#PLINK restrictions
#r2: LD analysis
#ld window kb: max distance of SNPs
#ld window r2: minimum value to report

#this PLINK command outputs a file with the LD results in a .ld extension

#fourth step: visualization of LD table for the complete binary data set
ld_analysis_complete <- read.table("LD_results/file_name.ld", header = TRUE)
write.csv(ld_analysis_complete, "LD_results/file_name.csv", row.names = FALSE)

#the .ld table column headers explain the values of each row
#the .CSV file has the row.names = FALSE to avoid creating a column with row
#numbers

#------------Complete LD analysis filtered by High and Low LD values -----------

#obtain high linked SNP pairs: filtered by R2 >= 0.3 & R2 <= 0.5
high_linked_complete <- ld_analysis_complete %>% filter(R2 >= 0.3 & R2 <= 0.5)

#obtain low linked SNP pairs: filtered by r2 >= 0.05
low_linked_complete <- ld_analysis_complete %>% filter(R2 >= 0.05 & R2 < 0.3)

#R2 it's the column name for the r2 value of LD
#you can put any values of r2 to filter in the function filter()

#save complete filtered LD analysis tables
write.csv(high_linked_complete, "LD_results/file_name.csv",row.names = FALSE) 
write.csv(low_linked_complete, "LD_results/file_name.csv", row.names = FALSE)

#----------------------------LD analysis filtered by sex -----------------------

#LD analysis for males only (sex code 1 in .fam files) with the pruned binary data set
system(paste0("./plink.exe --bfile file_name_input_prefix",
              " --filter-males --r2 --ld-window-kb 1000 --ld-window-r2 0.05",
              " --out LD_results/file_name_output_prefix"))

#LD analysis for females only (sex code 2 in .fam files) with the pruned binary data set
system(paste0("./plink.exe --bfile file_name_input_prefix",
              " --filter-females --r2 --ld-window-kb 1000 --ld-window-r2 0.05",
              " --out LD_results/file_name_output_prefix"))

#both PLINK commands output .ld files but filtered by sex

#-------------- Males LD analysis and filtration for High and Low LD--------------

ld_analysis_males <- read.table("LD_results/file_name.ld", header = TRUE)
write.csv(ld_analysis_males, "LD_results/file_name.csv", row.names = FALSE)

#obtain high linked SNP pairs: filtered by R2 >= 0.3 & R2 <= 0.5
high_linked_males <- ld_analysis_males %>% filter(R2 >= 0.3 & R2 <= 0.5)

#obtain low linked SNP pairs: filtered by R2 >= 0.05 & R2 < 0.3
low_linked_males <- ld_analysis_males %>% filter(R2 >= 0.05 & R2 < 0.3)

#save males filtered LD tables
write.csv(high_linked_males, "LD_results/file_name.csv",row.names = FALSE)
write.csv(low_linked_males, "LD_results/file_name.csv", row.names = FALSE)

#---------- Females LD analysis and filtration for High and Low LD ------------------

ld_analysis_females <- read.table("LD_results/file_name.ld", header = TRUE)
write.csv(ld_analysis_females, "LD_results/file_name.csv", row.names = FALSE)

#obtain high linked SNP pairs: filtered by R2 >= 0.3 & R2 <= 0.5
high_linked_females <- ld_analysis_females %>% filter(R2 >= 0.3 & R2 <= 0.5)

#obtain low linked SNP pairs: filtered by R2 >= 0.05 & R2 < 0.3
low_linked_females <- ld_analysis_females %>% filter(R2 >= 0.05 & R2 < 0.3)

#save females filtered LD tables for GWAS
write.csv(high_linked_females, "LD_results/file_name.csv", row.names = FALSE)
write.csv(low_linked_females, "LD_results/file_name.csv", row.names = FALSE)

# -------Table of SNPs with allele variants -------------------------------------

#this section of the script it's for creating a summary table with all the allelic
#variants of the SNPs that passed the pairwise filter

# read pruned .bim file to get allele information for the SNPS
# note: use the same prefix as your pruned dataset from step 2
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

# save summary table as csv
write.csv(snps_summary, "LD_results/file_name.csv", row.names = FALSE)

#print summary statistics
print(paste("Number of rows in SNP summary:", nrow(snps_summary)))
print(paste("Should equal pruned SNPs:", 3400))

#note: the expected number should match the number of SNPs in your pruned dataset
#if all the SNPs are biallelic. If not, modify the script to your necessities.