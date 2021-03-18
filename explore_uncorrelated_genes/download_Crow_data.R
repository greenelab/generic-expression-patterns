# sessionInfo
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base
#
#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.6        rstudioapi_0.13   magrittr_2.0.1    usethis_2.0.1     devtools_2.3.2    pkgload_1.1.0     R6_2.5.0
# [8] rlang_0.4.10      plyr_1.8.6        tools_3.6.3       pkgbuild_1.2.0    sessioninfo_1.1.1 cli_2.3.1         withr_2.4.1
# [15] ellipsis_0.3.1    remotes_2.2.0     assertthat_0.2.1  digest_0.6.27     rprojroot_2.0.2   lifecycle_1.0.0   crayon_1.4.1
# [22] processx_3.4.5    purrr_0.3.4       callr_3.5.1       fs_1.5.0          ps_1.5.0          curl_4.3          testthat_3.0.1
# [29] memoise_2.0.0     glue_1.4.2        compiler_3.6.3    desc_1.2.0        prettyunits_1.1.1

# This script is downloading data associated with Crow et. al. paper
# Following instructions from: https://github.com/PavlidisLab/gemmaAPI.R
install.packages('plyr')
# sudo apt-get install libcurl4-openssl-dev r-base
# sudo apt-get install libxml2-dev r-base
install.packages('covr')
install.packages('usethis')
install.packages('devtools')
devtools::install_github('PavlidisLab/gemmaAPI.R')

# Read in supplementary table containing experiment ids to download
supplementary_filename <- "pnas.1802973116.sd01.csv"
metadata <- as.matrix(read.csv(supplementary_filename, header=TRUE, row.names=1))

gse_ids <- metadata[,"External.ID"]

#length(gse_ids)
for (i in 1:3){
  print(gse_ids[i])
  # Download expression data for selected experiment ids
  data <- t(as.data.frame(gemmaAPI::datasetInfo(gse_ids[i], request='data', filter=FALSE)))

  colnames(data) <- data[1,]

  # Combine expression data per experiment id
  if (i == 1){
    all_data <- data
  }
  else{
    all_data <- plyr::rbind.fill(all_data, data)
  }

}

# Save
write.table(all_data, "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data.tsv", sep="\t", row.names=TRUE, col.names=TRUE)

