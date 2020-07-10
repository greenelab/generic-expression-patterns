"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts provide supporting functions to run analysis notebooks.
These scripts include: replacing ensembl gene ids with hgnc symbols. 
"""

import pandas as pd
import os


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


def subset_samples(samples_to_remove, num_runs, local_dir, project_id):
    """
    Removes user selected samples from the simulated experiments. This function
    overwrites the data in the simulated gene expression data files. 
    
    Arguments
    ---------
    samples_to_remove: lst
        list of samples ids to remove from each simulated experiment
    num_runs: int
        Number of simulated experiments
    local_dir: str
        Local directory containing simulated experiments
    project_id: str
        Project id to use to retrieve simulated experiments

    """

    for i in range(num_runs):
        simulated_data_file = os.path.join(
            local_dir,
            "pseudo_experiment",
            "selected_simulated_data_" + project_id + "_" + str(i) + ".txt",
        )

        # Read simulated data
        simulated_data = pd.read_csv(
            simulated_data_file, header=0, sep="\t", index_col=0
        )

        # Drop samples
        simulated_data = simulated_data.drop(samples_to_remove)

        # Save
        simulated_data.to_csv(simulated_data_file, float_format="%.5f", sep="\t")


def concat_simulated_data(local_dir, num_runs, project_id):
    """
    This function will concatenate the simulated experiments into a single dataframe
    in order to aggregate statistics across all experiments.

    Arguments
    ---------
    local_dir: str
        Local directory containing simulated experiments
    num_runs: int
        Number of simulated experiments
    project_id: str
        Project id to use to retrieve simulated experiments

    Returns
    -------
    Dataframe containing all simulated experiments concatenated together

    """

    simulated_DE_stats_all = pd.DataFrame()
    for i in range(num_runs):
        simulated_DE_stats_file = os.path.join(
            local_dir,
            "DE_stats",
            "DE_stats_simulated_data_" + project_id + "_" + str(i) + ".txt",
        )

        # Read results
        simulated_DE_stats = pd.read_csv(
            simulated_DE_stats_file, header=0, sep="\t", index_col=0
        )

        simulated_DE_stats.reset_index(inplace=True)

        # Concatenate df
        simulated_DE_stats_all = pd.concat([simulated_DE_stats_all, simulated_DE_stats])

    return simulated_DE_stats_all


def generate_summary_table(
    template_DE_stats, simulated_DE_summary_stats, col_to_rank, local_dir
):
    """
    Generate a summary table of the template and summary statistics

    Arguments
    ---------
    template_DE_stats: df
        dataframe containing DE statistics for template experiment
    simulated_DE_summary_stats: df
        dataframe containing aggregated DE statistics across all simulated experiments
    col_to_rank: str
        DE statistic to use to rank genes
    local_dir: str
        path to local machine where output file will be stored

    Returns
    -------
    Dataframe summarizing gene ranking for template and simulated experiments

    """
    # Merge template statistics with simulated statistics
    template_simulated_DE_stats = template_DE_stats.merge(
        simulated_DE_summary_stats, left_index=True, right_index=True
    )
    print(template_simulated_DE_stats.shape)

    # Parse columns
    median_pval_simulated = template_simulated_DE_stats[("adj.P.Val", "median")]
    mean_test_simulated = template_simulated_DE_stats[(col_to_rank, "mean")]
    std_test_simulated = template_simulated_DE_stats[(col_to_rank, "std")]
    count_simulated = template_simulated_DE_stats[(col_to_rank, "count")]
    rank_simulated = template_simulated_DE_stats[("ranking", "")]

    summary = pd.DataFrame(
        data={
            "Gene ID": template_simulated_DE_stats.index,
            "Adj P-value (Real)": template_simulated_DE_stats["adj.P.Val"],
            "Rank (Real)": template_simulated_DE_stats["ranking"],
            "Test statistic (Real)": template_simulated_DE_stats[col_to_rank],
            "Median adj p-value (simulated)": median_pval_simulated,
            "Rank (simulated)": rank_simulated,
            "Mean test statistic (simulated)": mean_test_simulated,
            "Std deviation (simulated)": std_test_simulated,
            "Number of experiments (simulated)": count_simulated,
        }
    )
    summary["Z score"] = (
        summary["Test statistic (Real)"] - summary["Mean test statistic (simulated)"]
    ) / summary["Std deviation (simulated)"]

    # Save file
    summary_file = os.path.join(local_dir, "gene_summary_table_" + col_to_rank + ".tsv")

    summary.to_csv(summary_file, float_format="%.5f", sep="\t")

    return summary


def merge_ranks_to_compare(
    your_summary_gene_ranks_df,
    reference_gene_ranks_file,
    reference_gene_name_col,
    reference_rank_col,
):
    """
    Given dataframes of your ranking of genes and reference ranking
    of genes. This function merges the ranking into one dataframe
    to be able to compare, `shared_gene_rank_df`

    Arguments
    ---------
    your_summary_gene_ranks_df: df
        dataframe containing your rank per gene
    reference_gene_ranks_file: file
        file contining reference ranks per gene
    reference_gene_name_col: str
        column header containing the reference genes
    reference_rank_col: str
        column header containing the reference rank

    Returns
    -------
    Dataframe containing your ranking and the reference ranking per gene

    """
    # Read in reference gene ranks file
    reference_gene_ranks_df = pd.read_csv(reference_gene_ranks_file, header=0, sep="\t")

    # Get list of our genes
    gene_ids = list(your_summary_gene_ranks_df.index)

    # Get list of published generic genes
    published_generic_genes = list(reference_gene_ranks_df[reference_gene_name_col])

    # Get intersection of gene lists
    shared_genes = set(gene_ids).intersection(published_generic_genes)

    # Get rank of shared genes
    your_gene_rank_df = pd.DataFrame(
        your_summary_gene_ranks_df.loc[shared_genes, "Rank (simulated)"]
    )

    # Merge published ranking
    shared_gene_rank_df = pd.merge(
        your_gene_rank_df,
        reference_gene_ranks_df[[reference_rank_col, reference_gene_name_col]],
        left_index=True,
        right_on=reference_gene_name_col,
    )

    return shared_gene_rank_df


def scale_reference_ranking(merged_gene_ranks_df, reference_rank_col):
    """
    In the case where the reference ranking and your ranking are not
    in the same range, this function scales the reference ranking
    to be in the range as your ranking.

    For example, if reference ranking ranged from (0,1) and your
    ranking ranged from (0,100). This function would scale the
    reference ranking to also be between 0 and 100.

    Note: This function is assuming that the reference ranking range
    is smaller than yours

    Arguments
    ---------
    merged_gene_ranks: df
        dataframe containing your rank and reference rank per gene
    reference_rank_col: str
        column header containing the reference rank

    Returns
    -------
    The same merged_gene_ranks dataframe with reference ranks re-scaled
    """
    # Scale published ranking to our range
    max_rank = max(merged_gene_ranks_df["Rank (simulated)"])
    merged_gene_ranks_df[reference_rank_col] = round(
        merged_gene_ranks_df[reference_rank_col] * max_rank
    )

    return merged_gene_ranks_df

