"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts perform calculations used in analysis notebooks.
These calculations include: getting spearman correlation and CI. 
"""

from scipy import stats


def spearman_ci(ci, gene_rank_df, num_permutations):
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
    """

    r, p = stats.spearmanr(
        gene_rank_df["Rank (simulated)"], gene_rank_df["DE_Prior_Rank"]
    )

    r_perm_values = []
    for i in range(num_permutations):

        sample = gene_rank_df.sample(n=len(gene_rank_df), replace=True)

        r_perm, p_perm = stats.spearmanr(
            sample["Rank (simulated)"], sample["DE_Prior_Rank"]
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


def aggregate_stats(col_to_rank, simulated_DE_stats_all):
    """
    Aggregate statistics across all simulated experiments

    Arguments
    ---------
    col_to_rank: str
        The DE statistic to use to rank genes. These are column headers of the DE
        statistics results table.
    simulated_DE_stats_all: df
        Dataframe of concatenated simulated experiments

    Returns
    --------
    Dataframe grouped by gene (row of the dataframe).
    For each gene this function calculates the median, mean, standard deviation
    of the distribution of the selected statistic (`col_to_rank`) across all
    simulated experiments.
    For each gene, it also returns the count to tell you the number of simulated
     experiments that were generated.
    """
    if col_to_rank == "adj.P.Val":
        simulated_DE_summary_stats = simulated_DE_stats_all.groupby(["index"])[
            [col_to_rank]
        ].agg(["median", "mean", "std", "count"])
    else:
        simulated_DE_summary_stats = simulated_DE_stats_all.groupby(["index"])[
            [col_to_rank, "adj.P.Val"]
        ].agg(
            {col_to_rank: ["median", "mean", "std", "count"], "adj.P.Val": ["median"]}
        )
    return simulated_DE_summary_stats


def rank_genes(col_to_rank, DE_summary_stats, is_template):
    """
    Returns the input dataframe (`DE_summary_stats`), ranked by the selected
    statistic, `col_to_rank` (if the input is the template experiment)
     or the median of the selected statistic (if the input is the simulated experiments).
    The ordering of the ranking depends on the statistic selected.

    """
    # If ranking by p-value or adjusted p-value then high rank = low value
    if col_to_rank in ["P.Value", "adj.P.Val"]:
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
    elif col_to_rank in ["logFC", "t"]:
        if is_template:
            DE_summary_stats["ranking"] = (
                DE_summary_stats[col_to_rank].abs().rank(ascending=True)
            )
            DE_summary_stats = DE_summary_stats.sort_values(
                by=col_to_rank, ascending=False
            )
        else:
            DE_summary_stats["ranking"] = (
                DE_summary_stats[(col_to_rank, "median")].abs().rank(ascending=True)
            )
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
