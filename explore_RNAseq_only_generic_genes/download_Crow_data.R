# sessionInfo
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#
#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base
#
# loaded via a namespace (and not attached):
#  [1] rstudioapi_0.13     magrittr_2.0.1      hms_1.0.0           tidyselect_1.1.0    R6_2.5.0            rlang_0.4.10
# [7] fansi_0.4.2         blob_1.2.1          httr_1.4.2          dplyr_1.0.5         tools_3.6.3         R.oo_1.24.0
# [13] utf8_1.2.1          DBI_1.1.0           gemmaAPI_2.0.3.9000 cli_2.3.1           ellipsis_0.3.1      assertthat_0.2.1
# [19] tibble_3.1.0        lifecycle_1.0.0     crayon_1.4.1        purrr_0.3.4         readr_1.4.0         vctrs_0.3.6
# [25] R.utils_2.10.1      curl_4.3            glue_1.4.2          compiler_3.6.3      pillar_1.5.1        generics_0.1.0
# [31] R.methodsS3_1.8.1   pkgconfig_2.0.3

# This script is downloading data associated with Crow et. al. paper
# Following instructions from: https://github.com/PavlidisLab/gemmaAPI.R
install.packages('plyr')
install.packages('knitr')
# sudo apt-get install libcurl4-openssl-dev r-base
# sudo apt-get install libxml2-dev r-base
install.packages('covr')
install.packages('usethis')
install.packages('devtools')
devtools::install_github('PavlidisLab/gemmaAPI.R')

# Read in supplementary table containing experiment ids to download
# Note manually adjustments:
# Change name: GSE11512-Human --> GSE11512.1
# Remove "GSE25119"; "GSE33903"; "GSE36809", "GSE7214"
supplementary_filename <- "/home/alexandra/Documents/Repos/generic-expression-patterns/explore_uncorrelated_genes/pnas.1802973116.sd01.csv"
metadata <- as.matrix(read.csv(supplementary_filename, header=TRUE, row.names=1))

gse_ids <- metadata[,"External.ID"]

for (i in 1:length(gse_ids)){
  print(gse_ids[i])
  # Download expression data for selected experiment ids
  data <- as.data.frame(gemmaAPI::datasetInfo(gse_ids[i], request='data', filter=FALSE))

  # Drop all columns except gene symbol and sample ids
  columns_to_drop <- c("Probe", "Sequence", "GeneName", "GemmaId", "NCBIid")
  data_tmp1 <- data[, !colnames(data) %in% columns_to_drop]

  # Remove NA column headers
  data_tmp2 <- data_tmp1[!is.na(data_tmp1$GeneSymbol),]

  # Only include unique genes because I'm not sure how to aggregate across multiple gene measurements
  data_tmp3 <- subset(data_tmp2, !duplicated(data_tmp2$GeneSymbol))

  # Get duplicate gene ids since duplicated() leaves a single instance of the duplicated gene
  duplicate_gene_ids <- unique(data_tmp2[duplicated(data_tmp2$GeneSymbol), "GeneSymbol"])

  # Set index to be gene symbols
  rownames(data_tmp3) <- data_tmp3$GeneSymbol

  # Drop extra `GeneSymbol row`
  data_tmp3 <- data_tmp3[, !colnames(data_tmp3) == "GeneSymbol"]

  # Transpose matrix to be sample x gene so we care merge on columns
  data_tmp4 <- as.data.frame(t(data_tmp3))

  # Remove gene ids that were duplicates, since duplicated() leaves a single instance of the gene
  data_final <- data_tmp4[, !colnames(data_tmp4) %in% duplicate_gene_ids]

  # Combine expression data per experiment id
  if (i == 1){
    all_data <- data_final
  }
  else{
    #shared_gene_ids <- intersect(names(all_data), names(data_final))
    #all_data_intersect <- rbind(all_data[, shared_gene_ids], data_final[, shared_gene_ids])
    all_data <- dplyr::bind_rows(all_data, data_final)
  }
}

# Save
write.table(all_data, "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data_union.tsv", sep="\t")

# tmp call
#result_pt1_filename <- "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data_part1.tsv"
#result_pt2_filename <- "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data_part2.tsv"
#result_pt3_filename <- "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data_part3.tsv"

#result_part1 <- as.data.frame(read.csv(result_pt1_filename, sep="\t", header=TRUE, row.names=1))
#result_part2 <- as.data.frame(read.csv(result_pt2_filename, sep="\t", header=TRUE, row.names=1))
#result_part3 <- as.data.frame(read.csv(result_pt3_filename, sep="\t", header=TRUE, row.names=1))

#shared_gene_ids <- intersect(names(result_part1), names(result_part2))
#shared_gene_ids_final <- intersect(shared_gene_ids, names(result_part3))
#all_data_new <- rbind(result_part1[, shared_gene_ids_final], result_part2[, shared_gene_ids_final])
#all_data_new <- rbind(all_data_new[, shared_gene_ids_final], result_part3[, shared_gene_ids_final])

#write.table(all_data_new, "/home/alexandra/Documents/Data/Generic_expression_patterns/Crow_expression_data.tsv", sep="\t", row.names=TRUE, col.names=TRUE)
#----- Scratch work -------
# To confirm the logic of the above script runs as expected
#data1 <- as.data.frame(gemmaAPI::datasetInfo(gse_ids[1], request='data', filter=FALSE))
#data2 <- as.data.frame(gemmaAPI::datasetInfo(gse_ids[2], request='data', filter=FALSE))
#
# Drop all columns except gene symbol and sample ids
#columns_to_drop <- c("Probe", "Sequence", "GeneName", "GemmaId", "NCBIid")
#data1_tmp1 <- data1[, !colnames(data1) %in% columns_to_drop]
#data2_tmp1 <- data2[, !colnames(data2) %in% columns_to_drop]
#
# Remove NA column headers
#data1_tmp2 <- data1_tmp1[!is.na(data1_tmp1$GeneSymbol),]
#data2_tmp2 <- data2_tmp1[!is.na(data2_tmp1$GeneSymbol),]
#
# Only include unique genes because I'm not sure how to aggregate across multiple gene measurements
#data1_tmp3 <- subset(data1_tmp2, !duplicated(data1_tmp2$GeneSymbol))
#data2_tmp3 <- subset(data2_tmp2, !duplicated(data2_tmp2$GeneSymbol))
#
# Get duplicate gene ids since duplicated() leaves a single instance of the duplicated gene
#duplicate_gene_ids_data1 <- unique(data1_tmp2[duplicated(data1_tmp2$GeneSymbol), "GeneSymbol"])
#duplicate_gene_ids_data2 <- unique(data2_tmp2[duplicated(data2_tmp2$GeneSymbol), "GeneSymbol"])
#
# Set index to be gene symbols
#rownames(data1_tmp3) <- data1_tmp3$GeneSymbol
#rownames(data2_tmp3) <- data2_tmp3$GeneSymbol
#
# Drop extra `GeneSymbol row`
#data1_tmp3 <- data1_tmp3[, !colnames(data1_tmp3) == "GeneSymbol"]
#data2_tmp3 <- data2_tmp3[, !colnames(data2_tmp3) == "GeneSymbol"]
#
# Transpose matrix to be sample x gene so we care merge on columns
#data1_tmp4 <- as.data.frame(t(data1_tmp3))
#data2_tmp4 <- as.data.frame(t(data2_tmp3))
#
# Remove gene ids that were duplicates, since duplicated() leaves a single instance of the gene
#data1_tmp5 <- data1_tmp4[, !colnames(data1_tmp4) %in% duplicate_gene_ids_data1]
#data2_tmp5 <- data2_tmp4[, !colnames(data2_tmp4) %in% duplicate_gene_ids_data2]
#
# Combine with previous dataset
#shared_gene_ids <- intersect(names(data1_tmp5), names(data2_tmp5))
#new_data <- rbind(data1_tmp5[, shared_gene_ids], data2_tmp5[, shared_gene_ids])
#new_data <- merge(data1_tmp5, data2_tmp5, by=intersect(names(data1_tmp5), names(data2_tmp5)), all=TRUE)
#
# Check feature size is the same
#dim(new_data)[2] == length(intersect(names(data1_tmp5), names(data2_tmp5)))
#
# Check that number of samples sum
#dim(new_data)[1] == dim(data1_tmp5)[1] + dim(data2_tmp5)[1]