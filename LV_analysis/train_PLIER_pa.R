# Using R 3.6.3
# This script is run outside of the `generic_expression` conda environment since R3.6.3 is required
# to run PLIER and conda had issues installing this verion of R.

# sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
#
# Matrix products: default
#
# Random number generation:
# RNG:     Mersenne-Twister
# Normal:  Inversion
# Sample:  Rounding
#
#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#loaded via a namespace (and not attached):
#[1] compiler_3.6.3 tools_3.6.3

# Resources:
# https://github.com/greenelab/multi-plier/blob/master/util/plier_util.R
# https://greenelab.github.io/multi-plier/05-sle-wb_PLIER.nb.html#plier_model_training
# https://github.com/wgmao/PLIER/blob/master/R/Allfuncs.R

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