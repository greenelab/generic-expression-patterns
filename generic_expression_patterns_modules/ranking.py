"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts perform calculations used in analysis notebooks.
These calculations include: getting spearman correlation and CI. 
"""

from scipy import stats
import pandas as pd
from generic_expression_patterns_modules import process


def spearman_ci(ci, gene_rank_df, num_permutations, data_type):
    """
    Returns spearman correlation score and confidence interval
    
    Arguments
    ---------
    ci: float in range (0-1)
        A 95% or 0.95 confidence interval
    gene_ranking_df: df
        Dataframe containing the our rank and Crow et. al. rank
    num_permutations: int
        The number of permutations to estimate the confidence interval
    data_type: str
        Either 'DE' or 'GSEA'
    """
    if data_type.lower() == "de":
        r, p = stats.spearmanr(
            gene_rank_df["Rank (simulated)"], gene_rank_df["DE_Prior_Rank"]
        )
    elif data_type.lower() == "gsea":
        r, p = stats.spearmanr(
            gene_rank_df["Rank (simulated)"], gene_rank_df["Powers Rank"]
        )

    r_perm_values = []
    for i in range(num_permutations):

        sample = gene_rank_df.sample(n=len(gene_rank_df), replace=True)

        if data_type.lower() == "de":
            r_perm, p_perm = stats.spearmanr(
                sample["Rank (simulated)"], sample["DE_Prior_Rank"]
            )
        elif data_type.lower() == "gsea":
            r_perm, p_perm = stats.spearmanr(
                sample["Rank (simulated)"], sample["Powers Rank"]
            )
        r_perm_values.append(r_perm)

    alpha = 1 - ci

    sort_r_perm_values = sorted(r_perm_values)
    offset = int(num_permutations * (alpha / 2))

    return (
        r,
        p,
        sort_r_perm_values[offset],
        sort_r_perm_values[num_permutations - offset],
    )


def aggregate_stats(col_to_rank, simulated_stats_all, data_type):
    """
    Aggregate statistics across all simulated experiments

    Arguments
    ---------
    col_to_rank: str
        The DE statistic to use to rank genes. These are column headers of the DE
        statistics results table.
    simulated_stats_all: df
        Dataframe of concatenated simulated experiments
    data_type: str
        Either 'DE' or 'GSEA'

    Returns
    --------
    Dataframe grouped by gene (row of the dataframe).
    For each gene this function calculates the median, mean, standard deviation
    of the distribution of the selected statistic (`col_to_rank`) across all
    simulated experiments.
    For each gene, it also returns the count to tell you the number of simulated
    experiments that were generated.
    """
    if data_type.lower() == "gsea":
        if col_to_rank == "padj":
            simulated_summary_stats = simulated_stats_all.groupby(["pathway"])[
                [col_to_rank]
            ].agg(["median", "mean", "std", "count"])

        else:
            simulated_summary_stats = simulated_stats_all.groupby(["pathway"])[
                [col_to_rank, "padj"]
            ].agg({col_to_rank: ["median", "mean", "std", "count"], "padj": ["median"]})

    if data_type.lower() == "de":
        if "adj.P.Val" in simulated_stats_all.columns:
            if col_to_rank == "adj.P.Val":
                simulated_summary_stats = simulated_stats_all.groupby(["index"])[
                    [col_to_rank]
                ].agg(["median", "mean", "std", "count"])
            else:
                simulated_summary_stats = simulated_stats_all.groupby(["index"])[
                    [col_to_rank, "adj.P.Val"]
                ].agg(
                    {
                        col_to_rank: ["median", "mean", "std", "count"],
                        "adj.P.Val": ["median"],
                    }
                )
        else:
            if col_to_rank == "padj":
                simulated_summary_stats = simulated_stats_all.groupby(["index"])[
                    [col_to_rank]
                ].agg(["median", "mean", "std", "count"])
            else:
                simulated_summary_stats = simulated_stats_all.groupby(["index"])[
                    [col_to_rank, "padj"]
                ].agg(
                    {
                        col_to_rank: ["median", "mean", "std", "count"],
                        "padj": ["median"],
                    }
                )

    return simulated_summary_stats


def rank_genes_or_pathways(col_to_rank, DE_summary_stats, is_template):
    """
    Returns the input dataframe (`DE_summary_stats`) that has been modified such that
    genes are ranked by the selected statistic, `col_to_rank` 
    (if the input is the template experiment) or the median of the selected statistic
    (if the input is the simulated experiments).
    The ordering of the ranking depends on the statistic selected.

    Arguments
    ---------
    col_to_rank: str
        DE statistic to use to rank genes
    DE_summary_stats: df
        dataframe containing gene ranking for either template or simulated experiments
    is_template: bool
        if the DE_summary_stats df is for the template experiment or simulated experiments    

    """

    # If ranking by p-value or adjusted p-value then high rank = low value
    if col_to_rank in ["P.Value", "adj.P.Val", "pvalue", "padj"]:
        if is_template:
            DE_summary_stats["ranking"] = DE_summary_stats[col_to_rank].rank(
                ascending=False
            )
            DE_summary_stats = DE_summary_stats.sort_values(
                by=col_to_rank, ascending=True
            )
        else:
            DE_summary_stats["ranking"] = DE_summary_stats[
                (col_to_rank, "median")
            ].rank(ascending=False)
            DE_summary_stats = DE_summary_stats.sort_values(
                by=(col_to_rank, "median"), ascending=True
            )

    # If ranking by logFC then high rank = high abs(value)
    elif col_to_rank in ["logFC", "t", "log2FoldChange"]:
        if is_template:
            DE_summary_stats["ranking"] = DE_summary_stats[col_to_rank].rank(
                ascending=True
            )
            DE_summary_stats = DE_summary_stats.sort_values(
                by=col_to_rank, ascending=False
            )
        else:
            DE_summary_stats["ranking"] = DE_summary_stats[
                (col_to_rank, "median")
            ].rank(ascending=True)
            DE_summary_stats = DE_summary_stats.sort_values(
                by=(col_to_rank, "median"), ascending=False
            )

    # If ranking by Z-score then high rank = high value
    else:
        if is_template:
            DE_summary_stats["ranking"] = DE_summary_stats[col_to_rank].rank(
                ascending=True
            )
            DE_summary_stats = DE_summary_stats.sort_values(
                by=col_to_rank, ascending=False
            )
        else:
            DE_summary_stats["ranking"] = DE_summary_stats[
                (col_to_rank, "median")
            ].rank(ascending=True)
            DE_summary_stats = DE_summary_stats.sort_values(
                by=(col_to_rank, "median"), ascending=False
            )

    return DE_summary_stats


def process_and_rank_genes_pathways(
    template_stats_filename,
    local_dir,
    num_simulated_experiments,
    project_id,
    analysis_type,
    col_to_rank_by,
):
    """
    This function uses DE or GSEA statistics to rank genes or pathways
    by genericness (i.e. how changed a gene or pathway is in the
    template experiment or how changed a gene or pathway is across
    the set of simuulated experiments).

    For the template experiment,
    1. Take the absolute value for DE or GSEA statistics where we only care
    about the change and not the direction in terms of ranking
    (i.e. log2 fold change and t-statistic)
    2. Rank genes or pathways by the user defined statistic

    For simulated experiments,
    1. Aggregate DE or GSEA statistics across all experiments
    2. Take the absolute value for DE or GSEA statistics where we only care
    about the change and not the direction in terms of ranking
    (i.e. log2 fold change and t-statistic)
    3. Rank genes or pathways by the user defined statistic

    Arguments
    ----------
    template_stats_filename: str
        File containing DE or GSEA statistics
    local_dir: str
        path to local machine where output file will be stored
    num_simulated_experiments: int
        Number of experiments simulated
    project_id: str
        Experiment identifier
    analysis_type: 'DE' or 'GSEA'
    col_to_rank_by: str
        Statistic to use to rank genes or pathways by

    """
    # For template experiment

    # Read template experiment
    template_stats = pd.read_csv(
        template_stats_filename, sep="\t", index_col=0, header=0
    )

    # Take absolute value of logFC and t statistic
    template_stats = process.abs_value_stats(template_stats)

    # Rank genes in template experiment
    template_stats = rank_genes_or_pathways(col_to_rank_by, template_stats, True)

    # For simulated experiments

    # Concatenate simulated experiments
    simulated_stats_all = process.concat_simulated_data(
        local_dir, num_simulated_experiments, project_id, analysis_type
    )
    # Take absolute value of logFC and t statistic
    simulated_stats_all = process.abs_value_stats(simulated_stats_all)

    # Aggregate statistics across all simulated experiments
    simulated_summary_stats = aggregate_stats(
        col_to_rank_by, simulated_stats_all, analysis_type
    )
    # Rank genes in simulated experiments
    simulated_summary_stats = rank_genes_or_pathways(
        col_to_rank_by, simulated_summary_stats, False
    )

    return template_stats, simulated_summary_stats

