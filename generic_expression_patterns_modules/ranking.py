"""
Author: Alexandra Lee
Date Created: 16 June 2020

This script provide supporting functions to run analysis notebooks.

This script includes functions to format intermediate data files 
to prepare to compare gene/pathway ranking:

* function to calculate Spearman correlation
* function to concatenate simulated data results
* function to get absolute value of test statistics to use for ranking
* function to generate summary data files
* function to scale ranking

"""
import os
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler


def concat_simulated_data(local_dir, num_runs, project_id, data_type):
    """
    This function will concatenate the simulated experiments into a single dataframe
    in order to aggregate statistics across all simulated experiments.

    Arguments
    ---------
    local_dir: str
        Local directory containing simulated experiments
    num_runs: int
        Number of simulated experiments
    project_id: str
        Project id to use to retrieve simulated experiments
    data_type: 'DE' or 'GSEA'

    Returns
    -------
    Dataframe containing all simulated experiments concatenated together

    """

    simulated_stats_all = pd.DataFrame()
    for i in range(num_runs):
        if data_type.lower() == "de":
            simulated_stats_file = os.path.join(
                local_dir,
                "DE_stats",
                "DE_stats_simulated_data_" + project_id + "_" + str(i) + ".txt",
            )
        elif data_type.lower() == "gsea":
            simulated_stats_file = os.path.join(
                local_dir,
                "GSEA_stats",
                "GSEA_stats_simulated_data_" + project_id + "_" + str(i) + ".txt",
            )

        # Read results
        simulated_stats = pd.read_csv(
            simulated_stats_file, header=0, sep="\t", index_col=0
        )

        simulated_stats.reset_index(inplace=True)

        # Concatenate df
        simulated_stats_all = pd.concat([simulated_stats_all, simulated_stats])

    return simulated_stats_all


def abs_value_stats(simulated_DE_stats_all):
    """
    This function takes the absolute value of columns=[`logFC`, `t`, `NES`].
    For ranking genes or pathways, we only care about the magnitude of the change for
    the logFC, t, NES statistic, but not the direction.

    The ranking for each gene will be based on the mean absolute value of either
    logFC, t, NES statistic, depending on the user selection
    """
    if "logFC" in simulated_DE_stats_all.columns:
        simulated_DE_stats_all["logFC"] = simulated_DE_stats_all["logFC"].abs()
    elif "log2FoldChange" in simulated_DE_stats_all.columns:
        simulated_DE_stats_all["log2FoldChange"] = simulated_DE_stats_all[
            "log2FoldChange"
        ].abs()
    elif "t" in simulated_DE_stats_all.columns:
        simulated_DE_stats_all["t"] = simulated_DE_stats_all["t"].abs()
    elif "NES" in simulated_DE_stats_all.columns:
        simulated_DE_stats_all["NES"] = simulated_DE_stats_all["NES"].abs()
    return simulated_DE_stats_all


def generate_summary_table(
    template_stats_filename,
    template_ranking_summary,
    simulated_ranking_summary,
    col_to_rank,
    local_dir,
    pathway_or_gene,
    params,
):
    """
    Generate a summary table of that includes DE association or
    enrichment statistics for the template experiment and summary
    statistics for the simulated experiments. This table also
    includes gene or pathway ranking and a score to indicate
    specificity of each gene or pathway to the template in question.

    Arguments
    ---------
    template_stats_filename: str
        File containing DE or GSEA statistics for template experiment
    template_ranking_summary: df
        dataframe containing DE or GSEA statistics and ranking for template experiment
    simulated_ranking_summary: df
        dataframe containing aggregated DE or GSEA statistics across all simulated experiments
        and ranking
    col_to_rank: str
        DE statistic to use to rank genes
    local_dir: str
        path to local machine where output file will be stored
    pathway_or_gene: str
        'pathway' or 'gene' to indicate which analysis
    params: dict
        parameter dictionary read in by config file

    Returns
    -------
    Dataframe summarizing gene ranking for template and simulated experiments

    """
    # Read in template experiment
    template_stats = pd.read_csv(
        template_stats_filename, sep="\t", index_col=0, header=0
    )

    # Merge template statistics with simulated statistics
    template_simulated_summary_stats = template_ranking_summary.merge(
        simulated_ranking_summary, left_index=True, right_index=True
    )
    shared_genes = list(template_simulated_summary_stats.index)

    # Parse columns
    if "adj.P.Val" in template_simulated_summary_stats.columns:
        median_pval_simulated = template_simulated_summary_stats[
            ("adj.P.Val", "median")
        ]
        col_name = "adj.P.Val"
    else:
        median_pval_simulated = template_simulated_summary_stats[("padj", "median")]
        col_name = "padj"
    mean_test_simulated = template_simulated_summary_stats[(col_to_rank, "mean")]
    std_test_simulated = template_simulated_summary_stats[(col_to_rank, "std")]
    count_simulated = template_simulated_summary_stats[(col_to_rank, "count")]
    rank_simulated = template_simulated_summary_stats[("ranking", "")]

    # Get raw values for test statistic if we took the abs
    # for ranking.
    # If test statistic is either log2 fold change,
    # t-statistic, Normalized Enrichment Score
    # then we want to also report the raw
    # value for users to know the directionality
    abs_stats_terms = ["t", "NES", "logFC", "log2FoldChange"]

    # Set variable strings depends on analysis
    if pathway_or_gene.lower() == "pathway":
        index_header = "Pathway ID"
        test_statistic = params["rank_pathways_by"]
        if test_statistic in abs_stats_terms:
            test_statistic_label = f"abs({test_statistic})"

    elif pathway_or_gene.lower() == "gene":
        index_header = "Gene ID"
        test_statistic = params["rank_genes_by"]
        if test_statistic in abs_stats_terms:
            test_statistic_label = f"abs({test_statistic})"

    # Create summary table

    if test_statistic in abs_stats_terms:
        summary = pd.DataFrame(
            data={
                index_header: template_simulated_summary_stats.index,
                "Adj P-value (Real)": template_simulated_summary_stats[col_name],
                "Rank (Real)": template_simulated_summary_stats["ranking"],
                f"{test_statistic_label} (Real)": template_simulated_summary_stats[
                    col_to_rank
                ],
                f"{test_statistic} (Real)": template_stats.loc[
                    shared_genes, test_statistic
                ],
                "Median adj p-value (simulated)": median_pval_simulated,
                "Rank (simulated)": rank_simulated,
                f"Mean {test_statistic_label} (simulated)": mean_test_simulated,
                "Std deviation (simulated)": std_test_simulated,
                "Number of experiments (simulated)": count_simulated,
            }
        )
        summary["Z score"] = (
            summary[f"{test_statistic_label} (Real)"]
            - summary[f"Mean {test_statistic_label} (simulated)"]
        ) / summary["Std deviation (simulated)"]
    else:
        summary = pd.DataFrame(
            data={
                index_header: template_simulated_summary_stats.index,
                "Adj P-value (Real)": template_simulated_summary_stats[col_name],
                "Rank (Real)": template_simulated_summary_stats["ranking"],
                f"{test_statistic} (Real)": template_stats.loc[
                    shared_genes, test_statistic
                ],
                "Median adj p-value (simulated)": median_pval_simulated,
                "Rank (simulated)": rank_simulated,
                f"Mean {test_statistic} (simulated)": mean_test_simulated,
                "Std deviation (simulated)": std_test_simulated,
                "Number of experiments (simulated)": count_simulated,
            }
        )
        summary["Z score"] = (
            summary[f"{test_statistic} (Real)"]
            - summary[f"Mean {test_statistic} (simulated)"]
        ) / summary["Std deviation (simulated)"]

    return summary


def merge_ranks_to_compare(
    your_summary_ranks_df, reference_ranks_file, reference_name_col, reference_rank_col,
):
    """
    Given dataframes of simulation-based ranking of genes or pathways
    and reference ranking of genes or pathways.
    This function merges the ranking into one dataframe
    to be able to compare, `shared_gene_rank_df`

    Arguments
    ---------
    your_summary_ranks_df: df
        dataframe containing your rank per gene or pathway
    reference_ranks_file: file
        file contining reference ranks per gene or pathway
    reference_name_col: str
        column header containing the reference genes or pathway
    reference_rank_col: str
        column header containing the reference rank

    Returns
    -------
    Dataframe containing your ranking and the reference ranking per gene

    """
    # Read in reference ranks file
    reference_ranks_df = pd.read_csv(
        reference_ranks_file, header=0, index_col=0, sep="\t"
    )

    # Get list of our genes or pathways
    gene_or_pathway_ids = list(your_summary_ranks_df.index)

    # Get list of published generic genes or pathways
    if reference_name_col == "index":
        published_generic_genes_or_pathways = list(reference_ranks_df.index)
    else:
        published_generic_genes_or_pathways = list(
            reference_ranks_df[reference_name_col]
        )
    # Get intersection of gene or pathway lists
    shared_genes_or_pathways = set(gene_or_pathway_ids).intersection(
        published_generic_genes_or_pathways
    )

    # Get your rank of shared genes
    # Note: ranking was performed before intersection
    # So there will be some jumps in the ranking due to
    # genes that were not shared, but the ordering should
    # still be preserved and Spearman correlation will perform
    # its own ranking so there is no need to rerank
    your_rank_df = pd.DataFrame(
        your_summary_ranks_df.loc[shared_genes_or_pathways, "Rank (simulated)"]
    )

    # Merge published ranking
    if reference_name_col == "index":
        shared_gene_rank_df = pd.merge(
            your_rank_df,
            reference_ranks_df[[reference_rank_col]],
            left_index=True,
            right_index=True,
        )
    else:
        shared_gene_rank_df = pd.merge(
            your_rank_df,
            reference_ranks_df[[reference_rank_col, reference_name_col]],
            left_index=True,
            right_on=reference_name_col,
        )

    return shared_gene_rank_df


def scale_reference_ranking(merged_gene_ranks_df, reference_rank_col):
    """
    In the case where the reference ranking and simulation-based ranking are not
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
    scaler = MinMaxScaler(
        feature_range=(
            min(merged_gene_ranks_df["Rank (simulated)"]),
            max(merged_gene_ranks_df["Rank (simulated)"]),
        )
    )

    merged_gene_ranks_df[reference_rank_col] = scaler.fit_transform(
        np.array(merged_gene_ranks_df[reference_rank_col]).reshape(-1, 1)
    )

    return merged_gene_ranks_df


def get_shared_rank_scaled(
    summary_df, reference_filename, ref_gene_col, ref_rank_col, data_type
):
    """
    Returns `shared rank scaled` dataframe, which contains the ranking
    using the simulation-based approach and the reference manual method.
    This function also returns correlation values comparing those two rankings.

    Arguments
    ------------
    summary_df: dataframe
        Dataframe containing our ranking per gene along with other statistics
        associated with that gene.
    reference_filename: str
        File containing gene ranks from reference publication (Crow et. al.)
    ref_gene_col: str
        Name of column header containing reference gene symbols
    ref_rank_col: str
        Name of column header containing reference ranks of genes
    data_type: 'DE' or 'GSEA'

    Returns
    -------
    A tuple that includes two entries: the first is the shared rank scaled
    dataframe, the second is a dict of correlation values (r, p, ci_low, ci_high).

    """
    # Merge our ranking and reference ranking
    shared_rank_df = merge_ranks_to_compare(
        summary_df, reference_filename, ref_gene_col, ref_rank_col
    )

    if max(shared_rank_df["Rank (simulated)"]) != max(shared_rank_df[ref_rank_col]):
        shared_rank_scaled_df = scale_reference_ranking(shared_rank_df, ref_rank_col)
    else:
        shared_rank_scaled_df = shared_rank_df

    # Note: This is in case lowly expressed genes were not pre-filtered before DESeq
    # (Michael Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of
    # dispersion is difficult for genes with very, very low average counts. See the methods.
    # However, it doesn't really effect the outcome because these genes have almost no power for
    # detecting differential expression. Effects runtime though.
    shared_rank_scaled_df = shared_rank_scaled_df[
        ~shared_rank_scaled_df["Rank (simulated)"].isna()
    ]

    # Get correlation
    r, p, ci_low, ci_high = spearman_ci(0.95, shared_rank_scaled_df, 1000, data_type)

    correlations = {"r": r, "p": p, "ci_low": ci_low, "ci_high": ci_high}

    # Print out correlation values
    for k, v in correlations.items():
        print(k, "=", v)

    return (shared_rank_scaled_df, correlations)


def compare_gene_ranking(
    summary_df, reference_filename, ref_gene_col, ref_rank_col, output_figure_filename
):
    """
    Compare gene ranking and generate a SVG figure.
    Returns correlations to make debugging easier.

    Arguments
    ------------
    summary_df: dataframe
        Dataframe containing our ranking per gene along with other statistics
        associated with that gene.
    reference_filename: str
        File containing gene ranks from reference publication (Crow et. al.)
    ref_gene_col: str
        Name of column header containing reference gene symbols
    ref_rank_col: str
        Name of column header containing reference ranks of genes
    output_figure_filename: str
        Filename of output figure
    """

    shared_gene_rank_scaled_df, correlations = get_shared_rank_scaled(
        summary_df, reference_filename, ref_gene_col, ref_rank_col, data_type="DE"
    )

    fig = sns.jointplot(
        data=shared_gene_rank_scaled_df,
        x="Rank (simulated)",
        y=ref_rank_col,
        kind="hex",
        marginal_kws={"color": "white"},
    )

    fig.set_axis_labels(
        "SOPHIE", "DE prior (Crow et. al. 2019)", fontsize=14, fontname="Verdana"
    )

    fig.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )

    return correlations


def compare_pathway_ranking(summary_df, reference_filename, output_figure_filename):
    """
    Compare pathway ranking.
    Returns correlations to make debugging easier.

    Arguments
    ------------
    summary_df: dataframe
        Dataframe containing our ranking per pathway along with other statistics associated with that pathway
    reference_filename:
        File containing pathway ranks from reference publication (Powers et. al.)
    output_figure_filename: str
            Filename of output figure
    """

    # Column headers for generic pathways identified by Powers et. al.
    ref_gene_col = "index"
    ref_rank_col = "Powers Rank"

    shared_pathway_rank_scaled_df, correlations = get_shared_rank_scaled(
        summary_df, reference_filename, ref_gene_col, ref_rank_col, data_type="GSEA"
    )

    fig = sns.scatterplot(
        data=shared_pathway_rank_scaled_df,
        x="Rank (simulated)",
        y=ref_rank_col,
        color="#15527d",
        s=50,
        alpha=0.7,
        edgecolors="face",
    )

    fig.set_xlabel("SOPHIE", fontsize=14, fontname="Verdana")
    fig.set_ylabel("Powers et. al. 2018", fontsize=14, fontname="Verdana")

    fig.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )

    return correlations


def add_pseudomonas_gene_name_col(summary_gene_ranks, base_dir):
    """
    This function adds a column to the input dataframe
    that contains the gene name corresponding to the
    pseudomonas gene id

    Arguments
    ---------
    summary_gene_ranks: df
        Dataframe of ranks and other statistics per gene
    base_dir: str
        Path to repository directiory
    """

    # Gene number to gene name file
    gene_name_filename = os.path.join(
        base_dir,
        "pseudomonas_analysis",
        "data",
        "metadata",
        "Pseudomonas_aeruginosa_PAO1_107.csv",
    )

    # Read gene number to name mapping
    gene_name_mapping = pd.read_table(
        gene_name_filename, header=0, sep=",", index_col=0
    )

    gene_name_mapping = gene_name_mapping[["Locus Tag", "Name"]]

    gene_name_mapping.set_index("Locus Tag", inplace=True)

    # Format gene numbers to remove extraneous quotes
    gene_number = gene_name_mapping.index
    gene_name_mapping.index = gene_number.str.strip('"')

    gene_name_mapping.dropna(inplace=True)

    # Remove duplicate mapping
    # Not sure which mapping is correct in this case
    # PA4527 maps to pilC and still frameshift type 4
    # fimbrial biogenesis protein PilC (putative pseudogene)
    gene_name_mapping = gene_name_mapping[
        ~gene_name_mapping.index.duplicated(keep=False)
    ]

    # Add gene names
    summary_gene_ranks["Gene Name"] = summary_gene_ranks["Gene ID"].map(
        gene_name_mapping["Name"]
    )

    return summary_gene_ranks


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
    data_type: 'DE' or 'GSEA'
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
    Aggregate DE or GSEA statistics across all simulated experiments
    after simulated experiments have been concatenated together.

    Arguments
    ---------
    col_to_rank: str
        The DE statistic to use to rank genes. These are column headers of the DE
        statistics results table.
    simulated_stats_all: df
        Dataframe of concatenated simulated experiments
    data_type: 'DE' or 'GSEA'

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
    genes are ranked by the selected statistic, `col_to_rank` (if the input is the 
    template experiment) or the median of the selected statistic
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
    template_stats = abs_value_stats(template_stats)

    # Rank genes in template experiment
    template_stats = rank_genes_or_pathways(col_to_rank_by, template_stats, True)

    # For simulated experiments

    # Concatenate simulated experiments
    simulated_stats_all = concat_simulated_data(
        local_dir, num_simulated_experiments, project_id, analysis_type
    )
    # Take absolute value of logFC and t statistic
    simulated_stats_all = abs_value_stats(simulated_stats_all)

    # Aggregate statistics across all simulated experiments
    simulated_summary_stats = aggregate_stats(
        col_to_rank_by, simulated_stats_all, analysis_type
    )
    # Rank genes in simulated experiments
    simulated_summary_stats = rank_genes_or_pathways(
        col_to_rank_by, simulated_summary_stats, False
    )

    return template_stats, simulated_summary_stats
