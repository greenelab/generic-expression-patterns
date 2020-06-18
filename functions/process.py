"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts provide supporting functions to run analysis notebooks.
These scripts include: replacing ensembl gene ids with hgnc symbols. 
"""

import pandas as pd


def replace_ensembl_ids(expression_df, gene_id_mapping):
    """
    Replaces ensembl gene ids with hgnc symbols 
    
    Arguments
    ---------
    expression_df: df
        gene expression data matrix (sample x gene)
    gene_id_mapping: df
        Dataframe mapping ensembl ids (used in DE_stats_file) to hgnc symbols,
        used in Crow et. al.
    """
    # Read in data
    gene_expression = expression_df

    # Different cases of many-many mapping between gene ids
    # count = 0
    # symbols = []
    # for symbol, group in gene_id_mapping.groupby("ensembl_gene_id"):
    #    if group['hgnc_symbol'].nunique() > 1:
    #        print(group)
    #        count += 1
    #        symbols.append(symbol)
    # count

    # Case 1: Ensembl ids are paralogs (geneA, geneA_PAR_Y) and map to the
    # same hgnc symbol. Homologous sequences are paralogous
    # if they were separated by a gene duplication event: if a gene in an
    # organism is duplicated to occupy two different positions in the same
    # genome, then the two copies are paralogous

    # Remove paralogs

    gene_expression = gene_expression.iloc[
        :, ~gene_expression.columns.str.contains("PAR_Y")
    ]

    # Case 2: Same ensembl ids are mapped to different gene symbol twice (CCL3L1, CCL3L3)
    # ENSG00000187510.7  ENSG00000187510    C12orf74
    # ENSG00000187510.7  ENSG00000187510     PLEKHG7
    # Manually map based on what is found on ensembl site
    manual_mapping = {
        "ENSG00000187510.7": "PLEKHG7",
        "ENSG00000230417.11": "LINC00595",
        "ENSG00000255374.3": "TAS2R45",
        "ENSG00000276085.1": "CCL3L1",
    }

    gene_expression.rename(manual_mapping, axis="columns", inplace=True)

    # Case 3: Some rows are duplicates
    # ENSG00000223773.7	ENSG00000223773	CD99P1
    # ENSG00000124334.17	ENSG00000124334	IL9R

    # Keep first occurence of duplicated ensembl ids
    gene_id_mapping = gene_id_mapping.loc[
        ~gene_id_mapping.index.duplicated(keep="first")
    ]

    # Replace ensembl ids with gene symbol
    gene_expression.columns = gene_expression.columns.map(
        gene_id_mapping["hgnc_symbol"]
    )

    # Remove rows where we couldn't map ensembl id to gene symbol
    gene_expression = gene_expression.iloc[:, gene_expression.columns != ""]
    gene_expression = gene_expression.iloc[:, gene_expression.columns.notnull()]

    # Save
    return gene_expression

