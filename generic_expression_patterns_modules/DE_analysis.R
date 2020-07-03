## Run this once to setup environment
## Used R 3.6.3
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")

# library('limma')

get_DE_stats <- function(metadata_file, 
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
  expression_data <- t(as.matrix(read.table(expression_file, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.table(metadata_file, sep="\t", header=TRUE, row.names=1))
  
  # NOTE: It make sure the metadata is in the same order 
  # as the column names of the expression matrix.
  group <- interaction(metadata[,1])
  
  mm <- model.matrix(~0 + group)
  
  ## DEGs of simulated data
  # lmFit expects input array to have structure: gene x sample
  # lmFit fits a linear model using weighted least squares for each gene:
  fit <- lmFit(expression_data, mm)
  
  # Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
  # Samples are grouped based on experimental condition
  # The variability of gene expression is compared between these groups
  if (grepl("SRP000762_groups", metadata_file)) {
    contr <- makeContrasts(groupDEX_treated - groupcontrol, levels = colnames(coef(fit)))
  } else if (grepl("SRP057087_groups", metadata_file)) {
    contr <- makeContrasts(grouplesion - groupnot_lesion, levels = colnames(coef(fit)))
  } else if (grepl("SRP012656_groups", metadata_file)) {
    contr <- makeContrasts(groupTumor - groupNormal, levels = colnames(coef(fit)))
  }
  
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

create_volcano <- function(expression_file,
                           experiment_id,
                           pval,
                           local_dir) {

    # This functioni generates a volcano plot using the output from
    # the DE analysis script 'get_DE_stats' and output it to local_dir

    # Read in expression data
    res <- read.table(expression_file, header=TRUE)

    threshold <- 0.05

    # Make a basic volcano plot
    f <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'logFC',
                         y = 'adj.P.Val',
                         xlim = c(-2,2),
                         pCutoff = threshold,
                         FCcutoff = 1,
                         pointSize = 1.0,
                         labSize = 3.0,
                         xlab=bquote(~Log[2]~ 'fold change'),
                         ylab=bquote(-~Log[10]~ 'FDR adj p-value')
    )

    # Save
    out_file = paste(local_dir, "volcano_template_data_", experiment_id,".png", sep="")  
    ggsave(out_file, plot = f, dpi=300)
                           }
