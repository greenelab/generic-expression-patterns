# Using R 3.6.3

# Manually load operons.rda and geneinfo.rda into R environment
# from https://github.com/greenelab/ADAGEpath/blob/master/data/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("greenelab/ADAGEpath")

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
ADAGEpath::visualize_gene_network(
  selected_signatures = subset_signatures,
  model = eADAGEmodel,
  gene_color_value = annot
)

ADAGEpath::visualize_gene_network(
  selected_signatures = names(all_signatures),
  model = eADAGEmodel,
  gene_color_value = annot
  )