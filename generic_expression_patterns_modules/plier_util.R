# Imported from https://github.com/greenelab/multi-plier/blob/master/util/plier_util.R
PLIERNewData <- function(exprs.mat, seed = 12345) {
  # A wrapper function for applying PLIER to a data set. We use the following
  # genesets that come with PLIER: bloodCellMarkersIRISDMAP, svmMarkers,
  # and canonicalPathways. We set the k parameter for the PLIER model by
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

  # load PLIER pathway and cell type data
  data(bloodCellMarkersIRISDMAP)
  data(svmMarkers)
  data(canonicalPathways)

  # combine the pathway data from PLIER
  all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,
                                   canonicalPathways)
  print(all.paths)
  # what genes are common to the pathway data and the expression matrix
  cm.genes <- PLIER::commonRows(all.paths, exprs.mat)

  # row normalize
  exprs.norm <- PLIER::rowNorm(exprs.mat)

  # what should we set the minimum k parameter to in PLIER? estimate the number
  # of PC for the SVD decomposition
  set.k <- PLIER::num.pc(exprs.norm[cm.genes, ])

  # PLIER main function + return results
  plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], all.paths[cm.genes, ],
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)

  return(plier.res)

}