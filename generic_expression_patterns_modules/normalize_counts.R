# Load libraries
library("DESeq2")

MRnorm_expression <- function(expression_filename, metadata_filename, expression_out_filename, size_out_filename) {
	# Read in data
	# Note the expression data is transposed to gene x sample in order to run Limma
	expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
	metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))

	print("Checking sample ordering...")
	print(all.equal(colnames(expression_data), rownames(metadata)))
	#print(colnames(expression_data))
	dds <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)
	rm(expression_data)
	dds <- estimateSizeFactors(dds, type="ratio")
	size_factor <- sizeFactors(dds)

	normalized_expression <- counts(dds, normalized=TRUE)
	rm(dds)
	print(dim(normalized_expression))

	# Save file
	write.table(t(normalized_expression), file = expression_out_filename, row.names = T, sep = "\t", quote = F)
	rm(normalized_expression)

	write.table(size_factor, file = size_out_filename, row.names = T, sep = "\t", quote = F)
}