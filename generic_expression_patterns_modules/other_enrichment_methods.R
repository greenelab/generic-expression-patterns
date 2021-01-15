# Load libraries
library("fgsea")
library("limma")
library("GSVA")
library("DEFormats")
library("edgeR")
library("DESeq2")
library("clusterProfiler")

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
  
  # Format index
  gene_ids <- row.names(expression_data)
  pathway_ind <- ids2indices(pathway_DB_data, gene_ids)

  # Call ROAST
  enrich_pathways <- mroast(y, index=pathway_ind, design, contrast=ncol(design), nrot=10000, adjust.method="BH")                       

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

  # Call CAMERA
  enrich_pathways <- camera(y, pathway_DB_data, design, contrast=ncol(design), nrot=10000)                       

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
                          kcdf="Poisson",
                          parallel.sz=1,
                          verbose=TRUE
  )                       
  return(as.data.frame(enrich_pathways))
}

find_enriched_pathways_ORA <- function(expression_filename,
                                      metadata_filename,
                                      pathway_DB_filename
                                      ) {

  # This function performs DE analysis using DESeq.
  # Expression data in expression_file are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # expression_file: str
  #   File containing gene expression data
  #
  # metadata_file: str
  #   File containing mapping between sample id and group

  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  #print(head(expression_data))

  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)

  deseq_object <- DESeq(ddset)

  # Note parameter settings:
  # `independentFilter=False`: We have turned off the automatic filtering, which
  # filter filter out those tests from the procedure that have no, or little 
  # chance of showing significant evidence, without even looking at their test statistic.
  # Typically, this results in increased detection power at the same experiment-wide 
  # type I error, as measured in terms of the false discovery rate. 
  # cooksCutoff=True (default): Cook's distance as a diagnostic to tell if a single sample
  # has a count which has a disproportionate impact on the log fold change and p-values.
  # These genes are flagged with an NA in the pvalue and padj columns
  deseq_results <- results(deseq_object, independentFiltering=FALSE)

  deseq_results_df <-  as.data.frame(deseq_results)

  # Get DEGs
  threshold=0.05
  backgrd_genes <- row.names(deseq_results_df)
  degs <- deseq_results_df[deseq_results_df[,'padj']<threshold & abs(deseq_results_df[,'log2FoldChange'])>1,]
  degs_name <- row.names(degs)

  # Get over-represented pathways
  pathway_DB_data <- read.gmt(pathway_DB_filename)

  enrich_pathways <- enricher(degs_name,
                              universe=backgrd_genes,
                              pvalueCutoff=0.05,
                              pAdjustMethod="BH",
                              TERM2GENE=pathway_DB_data[, c("ont", "gene")]
  )                     
  return(as.data.frame(summary(enrich_pathways)))
}
