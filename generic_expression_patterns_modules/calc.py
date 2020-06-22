"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts perform calculations used in analysis notebooks.
These calculations include: getting spearman correlation and CI. 
"""

import pandas as pd
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

