"""
Author: Alexandra Lee
Date Created: 18 December 2020

This script provide supporting functions to run analysis notebooks.

This script includes functions to perform latent variable analysis
using features from either multiPLIER or eADAGE models.
"""

from glob import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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

    generic_threshold: int (0,100)
        Threshold to use to define generic genes. Based on
        Percentile (simulated) column
    """
    print(summary_data.shape)

    # Generic genes
    ls_generic_genes = list(
        (
            summary_data[summary_data["Percentile (simulated)"] >= generic_threshold]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of generic genes: {len(ls_generic_genes)}")

    # Other (non-generic) genes
    ls_other_genes = list(
        (
            summary_data[summary_data["Percentile (simulated)"] < generic_threshold]
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
    that were included in the multiPLIER or eADAGE analyses. We want to make
    sure that our gene lists obtained from SOPHIE vs multiPLIER or eADAGE
    are consistent. This will prevent indexing by a gene that doesn't
    exist and resulting in NA values downstream.

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)
    """
    model_genes = list(LV_matrix.index)

    processed_dict_genes = {}
    for gene_label, ls_genes in dict_genes.items():
        ls_genes_processed = list(set(model_genes).intersection(ls_genes))

        processed_dict_genes[gene_label] = ls_genes_processed

    return processed_dict_genes


def get_nonzero_LV_coverage(dict_genes, LV_matrix):
    """
    This function counts the number of LVs that each
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
        LV_series = (LV_matrix.loc[ls_genes] != 0).sum(axis=1)

        dict_nonzero_coverage[gene_label] = LV_series

    return dict_nonzero_coverage


def get_highweight_LV_coverage(dict_genes, LV_matrix, quantile=0.9):
    """
    This function count the number of LVs that each
    gene contributes a lot to (i.e. has a high negative or positive
    weight contribution).
    This function returns a dictionary [gene id]: number of LVs

    Note: Using the quantile means that each LV has the same number
    of high weight values. Also here we are using a quantile cutoff
    since our distribution is not normal (exponential PDF)

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)

    quantile: float(0,1)
        Quantile to use to threshold weights. Default set to 90th quantile.
    """
    thresholds_per_LV = LV_matrix.quantile(quantile)

    # Manually checked that genes selected as high weight
    # are above threshold using below print statements
    # print(thresholds_per_LV)
    # print(LV_matrix)
    # print(
    #    LV_matrix.loc[
    #        (LV_matrix.abs() > thresholds_per_LV)["Node2"].values, "Node2"
    #    ]
    # )
    dict_highweight_coverage = {}
    for gene_label, ls_genes in dict_genes.items():

        LV_series = (LV_matrix.abs() > thresholds_per_LV).sum(axis=1)[ls_genes]

        dict_highweight_coverage[gene_label] = LV_series

    return dict_highweight_coverage


def get_highweight_LV_coverage_pseudomonas(dict_genes, LV_matrix):
    """
    This function count the number of LVs that each
    gene contributes a lot to (i.e. has a high negative or positive
    weight contribution).

    The high weight genes are determined based on the eADAGE paper
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5532071/).
    Though the method is described in an earlier paper
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5700673/).
    In this paper genes are considered high weight if their weight
    is at least 2.5 standard deviations from the mean since weights
    are normally distributed.

    This function returns a dictionary [gene id]: number of LVs

    Arguments
    ---------
    dict_genes: dict
        Dictionary mapping gene ids to label="generic", "other"

    LV_matrix: df
        Dataframe containing contribution of gene to LV (gene x LV matrix)
    """
    eADAGE_std_cutoff = 2.5

    mean_per_LV = LV_matrix.mean()

    std_per_LV = LV_matrix.std() * eADAGE_std_cutoff

    upper_threshold = mean_per_LV + std_per_LV
    lower_threshold = mean_per_LV - std_per_LV

    # Manually checked that genes selected as high weight
    # are above threshold using below print statements
    # print(upper_threshold)
    # print(lower_threshold)
    # print(LV_matrix.head(10))
    # print((LV_matrix > upper_threshold).head(10)
    # print((LV_matrix > upper_threshold).loc["PA0007"].sum())
    # print((LV_matrix > upper_threshold).sum(axis=1).head(10))

    dict_highweight_coverage = {}
    for gene_label, ls_genes in dict_genes.items():
        HW_pos = (LV_matrix > upper_threshold).sum(axis=1)[ls_genes]
        HW_neg = (LV_matrix < lower_threshold).sum(axis=1)[ls_genes]
        LV_series = HW_pos.add(HW_neg)

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


def get_prop_highweight_generic_genes(dict_genes, LV_matrix, quantile=0.9):
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

    quantile: float(0,1)
        Quantile to use to threshold weights. Default set to 90th quantile.
    """

    prop_highweight_generic_dict = {}
    generic_gene_ids = dict_genes["generic"]

    thresholds_per_LV = LV_matrix.quantile(quantile)
    # print(thresholds_per_LV)
    num_highweight_genes = (LV_matrix.abs() > thresholds_per_LV).sum()[0]

    # Manually checks
    # Note: all LV have the same number of total high weight genes since
    # we used quantile here
    # print((LV_matrix.abs() > thresholds_per_LV).sum())
    # print(num_highweight_genes)

    for LV_id in LV_matrix.columns:
        # print(thresholds_per_LV[LV_id])
        highweight_genes_per_LV = list(
            LV_matrix[(LV_matrix.abs() > thresholds_per_LV)[LV_id] == True].index
        )
        # print(LV_matrix.abs()[LV_id])
        # print((LV_matrix.abs() > thresholds_per_LV)[LV_id])
        # print(highweight_genes_per_LV)
        # break

        num_highweight_generic_genes = len(
            set(generic_gene_ids).intersection(highweight_genes_per_LV)
        )
        prop_highweight_generic_genes = (
            num_highweight_generic_genes / num_highweight_genes
        )
        prop_highweight_generic_dict[LV_id] = prop_highweight_generic_genes

    return prop_highweight_generic_dict


def get_prop_highweight_generic_genes_pseudomonas(dict_genes, LV_matrix):
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
    """
    eADAGE_std_cutoff = 2.5

    prop_highweight_generic_dict = {}
    generic_gene_ids = dict_genes["generic"]

    mean_per_LV = LV_matrix.mean()

    std_per_LV = LV_matrix.std() * eADAGE_std_cutoff

    upper_threshold = mean_per_LV + std_per_LV
    lower_threshold = mean_per_LV - std_per_LV

    num_highweight_pos_genes = (LV_matrix > upper_threshold).sum()
    num_highweight_neg_genes = (LV_matrix < lower_threshold).sum()
    num_highweight_genes = num_highweight_pos_genes.add(num_highweight_neg_genes)
    # print((LV_matrix > upper_threshold).sum())
    # print((LV_matrix < lower_threshold).sum())
    # print(num_highweight_genes)

    for LV_id in LV_matrix.columns:
        # print(LV_matrix[LV_id])
        # print(upper_threshold[LV_id])
        # print(lower_threshold[LV_id])
        pos_highweight_genes_per_LV = list(
            LV_matrix[(LV_matrix > upper_threshold)[LV_id] == True].index
        )
        neg_highweight_genes_per_LV = list(
            LV_matrix[(LV_matrix < lower_threshold)[LV_id] == True].index
        )
        highweight_genes_per_LV = (
            pos_highweight_genes_per_LV + neg_highweight_genes_per_LV
        )
        # print(pos_highweight_genes_per_LV)
        # print(neg_highweight_genes_per_LV)
        # print(highweight_genes_per_LV)
        # print(num_highweight_genes[LV_id])

        num_highweight_generic_genes = len(
            set(generic_gene_ids).intersection(highweight_genes_per_LV)
        )
        prop_highweight_generic_genes = (
            num_highweight_generic_genes / num_highweight_genes[LV_id]
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

    Note: This is only used for multiPLIER model, where we have
    information of LV and pathways associations.

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


def plot_dist_weights(
    LV_id, LV_matrix, shared_genes, num_genes, gene_id_mapping, out_filename
):
    """
    This function creates a distribution of weights for selected
    `LV_id`. This allows us to explore the contribution of genes
    to this LV

    Arguments
    ----------
    LV_id: str
        identifier for LV
    LV_matrix: df
        gene x LV matrix with weight values
    shared_genes: list
        list of genes that are shared by the multiPLIER or eADAGE analysis
        (so they have LV weight information) and SOPHIE analysis (so they have
        generic label)
    num_genes: int
        Number of genes to display
    gene_id_mapping: df
        dataframe containing mapping between genes and "generic" or "other"
        label
    out_filename: str
        file to save plot to
    """
    # Get index name
    LV_matrix.index.rename("geneID", inplace=True)

    # Get gene with num_gene top weights
    top_genes = list(LV_matrix.loc[shared_genes, LV_id].abs().nlargest(num_genes).index)
    weight_df = LV_matrix.loc[top_genes].reset_index()
    print(weight_df[LV_id])

    # Add label for if generic or not
    gene_ids = list(weight_df["geneID"].values)

    weight_df["gene type"] = list(gene_id_mapping.loc[gene_ids, "gene type"].values)

    fig = sns.barplot(
        data=weight_df,
        x=LV_id,
        y="geneID",
        hue="gene type",
        hue_order=["generic", "other"],
        dodge=False,
        palette=["#81448e", "lightgrey"],
    )
    L = plt.legend()
    L.get_texts()[0].set_text("Common")
    L.get_texts()[1].set_text("Other")

    fig.set_xlabel("Weight", fontsize=14, fontname="Verdana")
    fig.set_ylabel("Gene symbol", fontsize=14, fontname="Verdana")
    fig.set_title(f"Weight distribution for {LV_id}", fontsize=14, fontname="Verdana")

    for label in fig.get_yticklabels():
        label.set_style("italic")

    plt.xlim(0, 7)

    fig.figure.savefig(
        out_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


def plot_dist_weights_pseudomonas(
    LV_id, LV_matrix, shared_genes, num_genes, gene_id_mapping, out_filename
):
    """
    This function creates a distribution of weights for selected
    `LV_id`. This allows us to explore the contribution of genes
    to this LV.

    Here we are looking at only those HW genes identified using
    2.5 standard deviation from the mean weight at the `LV_id`

    Arguments
    ----------
    LV_id: str
        identifier for LV
    LV_matrix: df
        gene x LV matrix with weight values
    shared_genes: list
        list of genes that are shared by the multiPLIER or eADAGE analysis
        (so they have LV weight information) and SOPHIE analysis (so they have
        generic label)
    num_genes: int
        Number of genes to display
    gene_id_mapping: df
        dataframe containing mapping between genes and "generic" or "other"
        label
    out_filename: str
        file to save plot to
    """
    # Get weight for LV_id
    LV_id_weight = LV_matrix[LV_id]

    # Calculate thresholds
    eADAGE_std_cutoff = 2.5
    mean_weight = LV_id_weight.mean()
    std_weight = LV_id_weight.std() * eADAGE_std_cutoff
    upper_threshold = mean_weight + std_weight
    lower_threshold = mean_weight - std_weight

    # Get high weight genes
    HW_pos_genes = list(LV_id_weight[(LV_id_weight > upper_threshold).values].index)
    HW_neg_genes = list(LV_id_weight[(LV_id_weight < lower_threshold).values].index)

    HW_genes = HW_pos_genes + HW_neg_genes

    # Sort HW genes by abs weight
    sorted_HW_genes = list(
        LV_id_weight[HW_genes].abs().sort_values(ascending=False).index
    )[0:num_genes]

    # Get gene with num_gene top weights
    LV_matrix.index.rename("geneID", inplace=True)
    weight_df = LV_matrix.loc[sorted_HW_genes, LV_id].reset_index()
    print(weight_df)
    # Add label for if generic or not
    gene_ids = list(weight_df["geneID"].values)

    weight_df["gene type"] = list(gene_id_mapping.loc[gene_ids, "gene type"].values)

    fig = sns.barplot(
        data=weight_df,
        x=LV_id,
        y="geneID",
        hue="gene type",
        hue_order=["generic", "other"],
        dodge=False,
        palette=["#81448e", "lightgrey"],
    )

    fig.set_xlabel("Weight", fontsize=14, fontname="Verdana")
    fig.set_ylabel("Gene", fontsize=14, fontname="Verdana")
    fig.set_title(f"Weight distribution for {LV_id}", fontsize=14, fontname="Verdana")

    fig.figure.savefig(
        out_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )
