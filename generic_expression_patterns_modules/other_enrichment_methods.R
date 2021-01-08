# Load libraries
library("fgsea")
library("limma")
library("GSVA")

find_enriched_pathways_ROAST <- function(DE_stats_file,
                                   pathway_DB_filename,
                                   statistic) {

  # Read in data
  DE_stats_data <- read.table(DE_stats_file,
                              sep = "\t",
                              header = TRUE,
                              row.names = NULL)
  # Sort genes by feature 1

  # feature 1: numeric vector
  if (statistic == 'logFC') {
    col_num = 2
  } else if (statistic == 'log2FoldChange') {
    col_num = 3
  } else if (statistic == 't') {
    col_num = 4
  } else if (statistic == 'p-value') {
    col_num = 5
  } else if (statistic == 'adj p-value' || statistic == 'pvalue') {
    col_num = 6
  } else if ( statistic == 'padj') {
    col_num = 7
  }
  rank_genes <- as.numeric(as.character(DE_stats_data[, col_num]))

  # feature 2: named vector of gene ids
  names(rank_genes) <- as.character(DE_stats_data[,1])

  # feature 3: decreasing order
  rank_genes <- sort(rank_genes, decreasing = TRUE)

  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  enrich_pathways <- roast(pathways = pathway_DB_data,
                           stats = rank_genes,
                           nperm = 10000)                         

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
