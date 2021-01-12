# Load libraries
library("fgsea")
library("limma")
library("GSVA")
library("DEFormats")
library("edgeR")

find_enriched_pathways_ROAST <- function(expression_filename,
                                          metadata_filename,
                                          pathway_DB_filename) {

  # Read data
  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))
  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  # Create DEGList based on counts
  dge = DGEList(expression_data, group=group)

  # Design matrix
  design <- model.matrix(~0 + group)

  # Estimate dispersions
  y <- estimateDisp(dge, design)

  # Call ROAST
  enrich_pathways <- mroast(y, pathway_DB_data, design, contrast=ncol(design), nrot=1000)                       

  return(as.data.frame(enrich_pathways))
}

find_enriched_pathways_CAMERA <- function(expression_filename,
                                          metadata_filename,
                                          pathway_DB_filename) {

  # Read data
  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))
  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  # Create DEGList based on counts
  dge = DGEList(expression_data, group=group)

  # Design matrix
  design <- model.matrix(~0 + group)

  # Estimate dispersions
  y <- estimateDisp(dge, design)

  # Call ROAST
  enrich_pathways <- camera(y, pathway_DB_data, design, contrast=ncol(design), nrot=1000)                       

  return(as.data.frame(enrich_pathways))
}

find_enriched_pathways_GSVA <- function(expression_filename,
                                   pathway_DB_filename) {

  # Read in expression data
  # Transpose to gene x sample
  expression_data <- t(read.table(expression_filename,
                                  sep = "\t",
                                  header = TRUE,
                                  row.names = NULL))
  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  enrich_pathways <- gsva(expression_data,
                          pathway_DB_data,
                          kcdf="Poisson"
  )                       

  return(as.data.frame(enrich_pathways))
}
