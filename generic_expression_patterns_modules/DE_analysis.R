# Load libraries
library("limma")
library("DESeq2")

get_DE_stats_limma <- function(metadata_file,
                               experiment_id,
                               expression_file,
                               data_type,
                               local_dir,
                               run) {

  # This function performs DE analysis using expression data in expression_file
  # where samples are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # metadata_file: str
  #   File containing mapping between sample id and group
  #
  # experiment_id: str
  #   Experiment id used to label saved output filee
  #
  # expression_file: str
  #   File containing gene expression data
  #   Expression data should be of the form sample x gene
  #
  # data_type: str
  #   Either 'template' or 'simulated' to label saved output file
  #
  # local_dir: str
  #   Directory to save output files to
  #
  # run: str
  #   Used as identifier for different simulated experiments

  # Read in data
  # Note the expression data is transposed to gene x sample in order to run Limma
  expression_data <- t(
    as.matrix(
      read.csv(expression_file, sep="\t", header=TRUE, row.names=1)
    )
  )
  metadata <- as.matrix(
    read.csv(metadata_file, sep="\t", header=TRUE, row.names=1)
  )

  # NOTE: It make sure the metadata is in the same order
  # as the column names of the expression matrix.
  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  ## DEGs of simulated data
  # lmFit expects input array to have structure: gene x sample
  # lmFit fits a linear model using weighted least squares for each gene:
  fit <- lmFit(expression_data, mm)

  # Comparisons between groups (log fold-changes) are obtained as contrasts of
  # these fitted linear models:
  # Samples are grouped based on experimental condition
  # The variability of gene expression is compared between these groups
  contr <- makeContrasts(group2 - group1, levels = colnames(coef(fit)))

  # Estimate contrast for each gene
  tmp <- contrasts.fit(fit, contr)

  # Empirical Bayes smoothing of standard errors (shrinks standard errors
  # that are much larger or smaller than those from other genes towards the average standard error)
  tmp <- eBayes(tmp)

  # Get significant DEGs
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  all_genes <-  as.data.frame(top.table)

  # Find all DEGs based on adjusted p-value cutoff
  threshold <- 0.001
  num_sign_DEGs <- all_genes[all_genes[,'adj.P.Val']<threshold & abs(all_genes[,'logFC'])>1,]

  # Save summary statistics of DEGs
  if (data_type == "template") {
    out_file = paste(local_dir, "DE_stats/DE_stats_template_data_", experiment_id,"_", run, ".txt", sep="")
  } else if (data_type == "simulated") {
    out_file = paste(local_dir, "DE_stats/DE_stats_simulated_data_", experiment_id,"_", run, ".txt", sep="")
  }
  write.table(all_genes, file = out_file, row.names = T, sep = "\t", quote = F)

  return(nrow(num_sign_DEGs))

}

get_DE_stats_DESeq <- function(metadata_file,
                               experiment_id,
                               expression_file,
                               data_type,
                               local_dir,
                               run) {

  # This function performs DE analysis using DESeq.
  # Expression data in expression_file are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # metadata_file: str
  #   File containing mapping between sample id and group
  #
  # experiment_id: str
  #   Experiment id used to label saved output filee
  #
  # expression_file: str
  #   File containing gene expression data
  #
  # data_type: str
  #   Either 'template' or 'simulated' to label saved output file
  #
  # local_dir: str
  #   Directory to save output files to
  #
  # run: str
  #   Used as identifier for different simulated experiments

  expression_data <- t(as.matrix(read.csv(expression_file, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_file, sep="\t", header=TRUE, row.names=1))

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  #print(head(expression_data))

  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)

  deseq_object <- DESeq(ddset)

  deseq_results <- results(deseq_object)

  deseq_results_df <-  as.data.frame(deseq_results)

  # Save summary statistics of DEGs
  if (data_type == "template") {
    out_file = paste(local_dir, "DE_stats/DE_stats_template_data_", experiment_id,"_", run, ".txt", sep="")
  } else if (data_type == "simulated") {
    out_file = paste(local_dir, "DE_stats/DE_stats_simulated_data_", experiment_id,"_", run, ".txt", sep="")
  }
  write.table(deseq_results_df, file = out_file, row.names = T, sep = "\t", quote = F)
}
