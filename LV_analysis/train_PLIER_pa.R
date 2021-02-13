# Using R 3.6.3
# This script is run outside of the `generic_expression` conda environment since R3.6.3 is required
# to run PLIER and conda had issues installing this verion of R.

# Install PLIER
install.packages("remotes")
remotes::install_github("wgmao/PLIER")

expression_filename <- "https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl"

# Read in expression data
# PLIER expects input expression data to be of the form gene x sample
expression_data <- as.matrix(
  read.table(
    expression_filename,
    sep = "\t",
    header = TRUE,
    row.names = 1
    )
)

# Location of PLIER training script
source("../generic_expression_patterns/plier_util.R")

# run PLIER model
plier.result <- PLIERNewData(expression_data)

# TO DO: Save locally
model.file <- file.path("C:/Users/alexj/Documents/UPenn/CGreene/Remote/generic_expression_patterns/Pa_compendium_PLIER_model.RDS")
saveRDS(plier.result, file = model.file)