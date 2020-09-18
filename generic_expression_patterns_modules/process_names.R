# Functions to proces gene names

convert_ensembl_to_symbol <- function(ensembl.genes) {
  require(biomaRt)
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  getBM(
    values = ensembl.genes,
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters = 'ensembl_gene_id',
    mart = ensembl,
    bmHeader = FALSE
  )
}

get_ensembl_symbol_mapping <- function(input_file, out_file) {
  # Read in data
  data <- read.csv(input_file, sep = "\t", header = FALSE, row.names = 1)

  original_gene_id <- t(data[1,])

  # Remove version from gene id
  parsed_gene_id <- gsub("\\..*","", original_gene_id)

  # Get mapping from ensembl - gene symbol
  gene_id_mapping <- convert_ensembl_to_symbol(as.character(parsed_gene_id))

  write.table(
    gene_id_mapping,
    file = out_file,
    row.names = T,
    sep = "\t",
    quote = F
  )
}
