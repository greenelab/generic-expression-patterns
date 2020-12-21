"""
Author: Alexandra Lee
Date Created: 18 December 2020

This script provide supporting functions to run analysis notebooks.

This script includes functions to perform multiPLIER analysis
"""

from glob import glob
import pandas as pd


def get_generic_specific_genes(summary_data, generic_threshold):
    """
    This function returns a dictionary of generic genes and other
    (non-generic) genes, based on the statistics contained within the
    summary dataframes

    Here genes are determined as generic based on their
    ranking across multiple simulated experiments (i.e. generic genes
    are those that are high ranked = genes were found to be consistently
    changed across multiple simulated experiments. All other genes
    are 'other'

    Arguments
    ---------
    summary_data: df
        Dataframe containing gene summary statistics

    generic_threshold: int
        Threshold to use to define generic genes
    """
    print(summary_data.shape)

    # Generic genes
    ls_generic_genes = list(
        (
            summary_data[summary_data["Rank (simulated)"] >= generic_threshold]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of generic genes: {len(ls_generic_genes)}")

    # Other (non-generic) genes
    ls_other_genes = list(
        (
            summary_data[summary_data["Rank (simulated)"] < generic_threshold]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of other genes: {len(ls_other_genes)}")

    # Create dictionary
    dict_genes = {
        "generic": ls_generic_genes,
        "other": ls_other_genes,
    }

    return dict_genes


def process_generic_specific_gene_lists(dict_genes, LV_matrix):
    """
    This function returns the dictionary of generic genes and specific genes
    that were included in the multiplier analysis. 

    This prevents indexing by a gene that doesn't exist and resulting in NA values

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)
    """
    multiplier_genes = list(LV_matrix.index)

    processed_dict_genes = {}
    for gene_label, ls_genes in dict_genes.items():
        ls_genes_processed = list(set(multiplier_genes).intersection(ls_genes))

        processed_dict_genes[gene_label] = ls_genes_processed

    return processed_dict_genes


def get_nonzero_LV_coverage(dict_genes, LV_matrix):
    """
    This function count the number of LVs that each
    gene is present in (i.e. has a nonzero contribution).
    This function returns a dictionary [gene id]: number of LVs

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)
    """
    dict_nonzero_coverage = {}
    for gene_label, ls_genes in dict_genes.items():
        LV_series = (LV_matrix.loc[ls_genes] > 0).sum(axis=1)

        dict_nonzero_coverage[gene_label] = LV_series

    return dict_nonzero_coverage


def get_highweight_LV_coverage(
    dict_genes, LV_matrix, if_normalized=False, quantile=0.9
):
    """
    This function count the number of LVs that each
    gene contributes a lot to (i.e. has a high weight contribution).
    This function returns a dictionary [gene id]: number of LVs

    Note: If we didn't normalize per LV so each LV has the same number
    of high weight values.

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)

    if_normalized: bool
        True if LV_matrix is normalized per LV

    quantile: float(0,1)
        Quantile to use to threshold weights. Default set to 90th quantile.
    """
    if if_normalized:
        threshold = 0.063
        dict_highweight_coverage = {}
        for gene_label, ls_genes in dict_genes.items():
            LV_series = (LV_matrix > threshold).sum(axis=1)[ls_genes]

            dict_highweight_coverage[gene_label] = LV_series
    else:
        thresholds_per_LV = LV_matrix.quantile(quantile)

        dict_highweight_coverage = {}
        for gene_label, ls_genes in dict_genes.items():
            LV_series = (LV_matrix > thresholds_per_LV).sum(axis=1)[ls_genes]

            dict_highweight_coverage[gene_label] = LV_series

    return dict_highweight_coverage


def assemble_coverage_df(dict_genes, nonzero_dict, highweight_dict):
    """
    This function assembles the coverage dfs into
    one df to be used for plotting

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    nonzero_dict: dict
        Dictionary mapping [gene type]: number of LVs present

    highweight_dict: dict
        Dictionary mapping [gene type]: number of LVs gene is highweight in

    """
    all_coverage = []
    for gene_label in dict_genes.keys():
        merged_df = pd.DataFrame(
            nonzero_dict[gene_label], columns=["nonzero LV coverage"]
        ).merge(
            pd.DataFrame(
                highweight_dict[gene_label], columns=["highweight LV coverage"]
            ),
            left_index=True,
            right_index=True,
        )
        merged_df["gene type"] = gene_label
        all_coverage.append(merged_df)

    all_coverage_df = pd.concat(all_coverage)

    return all_coverage_df


def get_prop_highweight_generic_genes(
    dict_genes, LV_matrix, if_normalized=False, quantile=0.9
):
    """
    This function returns a dictionary mapping 
    [LV id]: proportion of high weight generic genes

    Arguments
    ---------
    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)
    
    if_normalized: bool
        True if LV_matrix is normalized per LV

    quantile: float(0,1)
        Quantile to use to threshold weights. Default set to 90th quantile.
    """

    prop_highweight_generic_dict = {}
    generic_gene_ids = dict_genes["generic"]

    if if_normalized:
        for LV_id in LV_matrix.columns:
            threshold = 0.063
            num_highweight_genes = (LV_matrix > threshold).sum()[0]
            highweight_genes_per_LV = list(
                LV_matrix[(LV_matrix > quantile)[LV_id] == True].index
            )

            num_highweight_generic_genes = len(
                set(generic_gene_ids).intersection(highweight_genes_per_LV)
            )
            prop_highweight_generic_genes = (
                num_highweight_generic_genes / num_highweight_genes
            )
            prop_highweight_generic_dict[LV_id] = prop_highweight_generic_genes
    else:
        thresholds_per_LV = LV_matrix.quantile(quantile)
        num_highweight_genes = (LV_matrix > thresholds_per_LV).sum()[0]

        for LV_id in LV_matrix.columns:
            highweight_genes_per_LV = list(
                LV_matrix[(LV_matrix > thresholds_per_LV)[LV_id] == True].index
            )

            num_highweight_generic_genes = len(
                set(generic_gene_ids).intersection(highweight_genes_per_LV)
            )
            prop_highweight_generic_genes = (
                num_highweight_generic_genes / num_highweight_genes
            )
            prop_highweight_generic_dict[LV_id] = prop_highweight_generic_genes

    return prop_highweight_generic_dict


def create_LV_df(
    prop_highweight_generic_dict,
    multiplier_model_summary,
    proportion_generic,
    out_filename,
):
    """
    This function creates and saves dataframe that contains the metadata
    associated with the LV that is contributed most by generic genes

    Arguments
    ---------
    prop_highweight_generic_dict: dict
        Dictionary mapping LV_id: proportion of generic genes that are high weight
    
    multiplier_model_summary: df
        Dataframe containing summary statistics for which pathways LV are significantly associated

    proportion_generic: float
        Threshold for the proportion of high weight genes to be generic in a LV
    """
    generic_LV = []
    for k, v in prop_highweight_generic_dict.items():
        if v > proportion_generic:
            print(k, v)
            generic_LV.append(k)

    if len(generic_LV) > 0:
        LV_ids = [int(i.replace("LV", "")) for i in generic_LV]

        LV_df = multiplier_model_summary[
            multiplier_model_summary["LV index"].isin(LV_ids)
        ]

        LV_df.to_csv(out_filename, sep="\t")
    else:
        print("No LVs with high proportion of generic genes")
