# Alternative enrichment analyses to GSEA
# These are methods were implemented, expecting
# RNA-seq input. These will need to be adjusted if
# using microarray data

# Load libraries
library("fgsea")
library("limma")
library("GSVA")
library("DEFormats")
library("edgeR")
library("DESeq2")
library("clusterProfiler")

find_enriched_pathways_CAMERA <- function(expression_filename,
                                          metadata_filename,
                                          pathway_DB_filename,
                                          platform) {

  # ---------------------------------------------------------
  # CAMERA (Correlation Adjusted MEan RAnk gene set test)
  # is based on the idea of estimating the variance inflation
  # factor associated with inter-gene correlation, and
  # incorporating this into parametric or rank-based
  # test procedures. (available in limma)
  # * Competitive gene set test
  # * Performs the same rank-based test procedure as GSEA,
  # but also estimates the correlation between genes,
  # instead of treating genes as independent
  # * Recall GSEA: 1) Rank all genes using DE association statistics.
  # 2) An enrichment score (ES) is defined as the maximum distance
  # from the middle of the ranked list. Thus, the enrichment score
  # indicates whether the genes contained in a gene set are clustered
  # towards the beginning or the end of the ranked list
  # (indicating a correlation with change in expression).
  # 3) Estimate the statistical significance of the ES by a
  # phenotypic-based permutation (permute samples assigned to label)
  # test in order to produce a null distribution for the ES
  # (i.e. scores based on permuted phenotype)
  # * Appropriate for small and large fold changes
  # (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/)
  # ---------------------------------------------------------

  # Read data
  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))
  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  # Design matrix
  design <- model.matrix(~0 + group)

  if (platform == "rnaseq"){
    # Create DGEList based on counts
    dge = DGEList(expression_data, group=group)

    # Estimate dispersions
    y <- estimateDisp(dge, design)

    # Call CAMERA
    enrich_pathways <- camera(y, pathway_DB_data, design, contrast=ncol(design), nrot=10000)
  }
  # else if data used microarray
  else{
    # Call CAMERA
    enrich_pathways <- camera(expression_data, pathway_DB_data, design, contrast=ncol(design), nrot=10000)
  }


  return(as.data.frame(enrich_pathways))
}

find_enriched_pathways_GSVA <- function(expression_filename,
                                        pathway_DB_filename,
                                        platform){

  # ---------------------------------------------------------
  # GSVA(Gene Set Variation Analysis) calculates sample-wise
  # gene set enrichment scores as a function of genes inside
  # and outside the gene set. This method is well-suited for
  # assessing gene set variation across a dichotomous phenotype.
  # (bioconductor package GSVA)
  # * Competitive gene set test
  # * Estimates variation of gene set enrichment over the samples
  # independently of any class label
  # (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/)
  #
  # Note: Gaussian CDF is suitable for continuous expression
  # values in microarrays, and a Poisson CDF is better for
  # counts in RNA-seq data).
  # https://www.bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf
  # ---------------------------------------------------------

  # Read in expression data
  # Transpose to gene x sample
  expression_data <- t(read.table(expression_filename,
                                  sep = "\t",
                                  header = TRUE,
                                  row.names = NULL))
  pathway_DB_data <- gmtPathways(pathway_DB_filename)

  if (platform == "rnaseq"){
    enrich_pathways <- gsva(expression_data,
                          pathway_DB_data,
                          kcdf="Poisson",
                          parallel.sz=1,
                          verbose=TRUE
  )
  }
  # else if data used microarray
  else{
    enrich_pathways <- gsva(expression_data,
                          pathway_DB_data,
                          kcdf="Gaussian",
                          parallel.sz=1,
                          verbose=TRUE
  )
  }

  return(as.data.frame(enrich_pathways))
}

find_enriched_pathways_ORA <- function(expression_filename,
                                      metadata_filename,
                                      pathway_DB_filename,
                                      platform) {
  # ---------------------------------------------------------
  # ORA (over-representation analysis) uses the hypergeometric
  # test to determine if there a significant over-representation
  # of pathway in the selected set of DEGs. Here we're using
  # clusterProfiler library but there are multiple options for this analysis.
  # (https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enricher)
  # ---------------------------------------------------------

  # Read data
  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  if (platform == "rnaseq"){

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
  }
  # else if data used microarray
  else{
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
    backgrd_genes <- row.names(all_genes)

    # Find all DEGs based on adjusted p-value cutoff
    threshold <- 0.05
    degs <- all_genes[all_genes[,'adj.P.Val']<threshold & abs(all_genes[,'logFC'])>1,]
    degs_name <- row.names(degs)
  }

  # Get over-represented pathways
  pathway_DB_data <- read.gmt(pathway_DB_filename)

  if (nrow(degs)==0){
    enrich_pathways <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)
  }
  else {
    enrich_pathways <- enricher(degs_name,
                                universe=backgrd_genes,
                                pvalueCutoff=1.0,
                                pAdjustMethod="BH",
                                qvalueCutoff = 1.0,
                                TERM2GENE=pathway_DB_data[, c("ont", "gene")]
    )
  }
  return(as.data.frame(summary(enrich_pathways)))
}
