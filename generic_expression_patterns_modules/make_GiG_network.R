# Using R 3.6.3

# Manually load operons.rda and geneinfo.rda into R environment
# from https://github.com/greenelab/ADAGEpath/blob/master/data/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("impute")) install.packages("impute")
if (!require("greenelab/ADAGEpath")) install.packages("greenelab/ADAGEpath")

# Get eADAGE model
eADAGEmodel <- ADAGEpath::eADAGEmodel

# Get all eADAGE signatures
all_signatures <- ADAGEpath::extract_signatures(ADAGEpath::eADAGEmodel)
subset_signatures <- c("Node1pos", "Node1neg")

# Read in dataframe mapping gene ids to
# generic=1 or not generic=0
annot_filename <- "annot_df.tsv"
annot <- as.data.frame(
  read.table(
    annot_filename,
    sep="\t",
    header=TRUE,
    row.names=NULL
    )
)
names(annot) <- c("geneID", "other")

# Plot G-G network
# In this network nodes = genes and edges = similar weight profiles for
# how much that gene contributes to the eADAGE (denoising autoencoder) latent variable.
# This network is showing all 5,549 measured genes.
# The clustering is using pearson correlation with a default cutoff of 0.5.

# Code below is in case we decide to only show edges based on a subset
# of LVs instead of all of them.
#ADAGEpath::visualize_gene_network(
#  selected_signatures = subset_signatures,
#  model = eADAGEmodel,
#  gene_color_value = annot
#)

ADAGEpath::visualize_gene_network(
  selected_signatures = names(all_signatures),
  model = eADAGEmodel,
  gene_color_value = annot
  )