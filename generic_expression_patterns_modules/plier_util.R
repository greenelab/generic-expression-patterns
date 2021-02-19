# Imported from https://github.com/greenelab/multi-plier/blob/master/util/plier_util.R
PLIERNewData <- function(exprs.mat, seed = 12345) {
  # A function for applying PLIER to a exprs.mat. We use
  # P. aeruginosa genesets from KEGG.
  # We set the k parameter for the PLIER model by
  # identifying the number of "significant PCs" with PLIER::num.pc and then
  # using sig PCs * 0.3. This is consistent with recommendations from the
  # PLIER authors.
  #
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   seed: an integer to be supplied to set.seed() for reproducibility
  #         purposes, default is 12345
  #
  # Returns:
  #   plier.res: output from PLIER::PLIER()
  #
  require(PLIER)

  set.seed(seed)

  # load Pa pathway data generated from
  # 0_process_Pa_KEGG_pathways.ipynb

  pa_paths_filename <- "..LV_analysis/pa_kegg_pathway_processed.tsv"
  pa.paths <- as.matrix(
    read.table(
      pa_pathway_filename,
      sep = "\t",
      header = TRUE,
      row.names = 1
      )
    )

  # what genes are common to the pathway data and the expression matrix
  cm.genes <- PLIER::commonRows(pa.paths, exprs.mat)

  # row normalize
  exprs.norm <- PLIER::rowNorm(exprs.mat)

  # what should we set the minimum k parameter to in PLIER? estimate the number
  # of PC for the SVD decomposition
  set.k <- PLIER::num.pc(exprs.norm[cm.genes, ])

  # PLIER main function + return results
  plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], pa.paths[cm.genes, ],
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)

  return(plier.res)

}