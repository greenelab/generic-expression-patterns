## Run this once to setup environment
## Used R 3.6.3
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

#library(clusterProfiler)

find_enriched_pathways <- function(DE_stats_file,
                                   pathway_DB,
                                   statistic){
    # Read in data
    DE_stats_data <- read.table(DE_stats_file, sep="\t", header=TRUE, row.names=NULL)
   
    # Sort genes by feature 1
    
    # feature 1: numeric vector
	if (statistic == 'logFC'){
		col_num = 2
	} else if (statistic == 'log2FoldChange'){
		col_num = 3
    } else if (statistic =='t'){
		col_num = 4
    } else if (statistic == 'p-value'){
		col_num = 5
	} else if (statistic == 'adj p-value' || statistic == 'pvalue'){
		col_num = 6
    } else if ( statistic == 'padj'){
		col_num = 7
    }
    rank_genes <- as.numeric(as.character(DE_stats_data[,col_num]))

    # feature 2: named vector of gene ids
    names(rank_genes) <- as.character(DE_stats_data[,1])

	## feature 3: decreasing order
	rank_genes <- sort(rank_genes, decreasing = TRUE)

    pathway_DB_data <- gmtPathways(hallmark_DB_file)
 
    #enrich_pathways <- GSEA(geneList=rank_genes, 
    #                        TERM2GENE=pathway_DB_data,
    #                        nPerm=100000,
    #                        by='fgsea',
    #                        verbose=T)
    enrich_pathways <- fgsea(pathways=pathway_DB_data,
                              stats=rank_genes,
                              nperm=10000)

    return(as.data.frame(enrich_pathways))
}