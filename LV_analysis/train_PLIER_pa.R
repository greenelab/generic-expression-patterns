# Using R 3.6.3
# This script was run outside of the `generic_expression` conda environment since R3.6.3 is required
# and conda had issues installing it

# Install PLIER
install.packages("remotes")
remotes::install_github("wgmao/PLIER")

# Location of PLIER training script
source("/home/alexandra/Documents/Repos/generic-expression-patterns/generic_expression_patterns_modules/plier_util.R")

expression_filename <- "https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl"
# Read in expression data
# Transpose to gene x sample
expression_data <- read.table(expression_filename,
                                sep = "\t",
                                header = TRUE,
                                row.names = NULL)

# run PLIER model
plier.result <- PLIERNewData(expression_data)

# Check
# Check expression input dimension and data type and note
# Check pathway annotation format
# Format KEGG pathway and load into plier_util.R
# Run and check output, particularly Z