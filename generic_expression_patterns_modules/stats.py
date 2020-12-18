"""
Author: Alexandra Lee
Date Created: 18 December 2020

This script provide supporting functions to run analysis notebooks.

This script includes functions to prepare the data to run
DE and GSEA analyses
"""
import os
import csv
import pandas as pd


def compare_and_reorder_samples(expression_file, metadata_file):
    """
    This function checks that the ordering of the samples matches
    between the expression file and the metadata file. This
    ordering is used for calculating DEGs.
    """
    # Check ordering of sample ids is consistent between gene expression data and metadata
    metadata = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)
    metadata_sample_ids = metadata.index

    expression_data = pd.read_csv(expression_file, sep="\t", header=0, index_col=0)
    expression_sample_ids = expression_data.index

    if metadata_sample_ids.equals(expression_sample_ids):
        print("sample ids are ordered correctly")
    else:
        # Convert gene expression ordering to be the same as
        # metadata sample ordering
        print("sample ids don't match, going to re-order gene expression samples")
        expression_data = expression_data.reindex(metadata_sample_ids)
        expression_data.to_csv(expression_file, sep="\t")


def format_pseudomonas_pathway_DB(pathway_DB_filename, local_dir, out_filename):
    """
    This function reads in pseudomonas pathway data from
    `pathway_DB_filename` and formats and outputs it
    to `output_filename` in order to be
    used in GSEA_analysis.R

    Note: Currently this function is specifically
    customized to expect pathway_DB_filename = 
    "https://raw.githubusercontent.com/greenelab/adage/master/Node_interpretation/pseudomonas_KEGG_terms.txt"

    """
    # Read in pathway data
    pa_pathway_DB = pd.read_csv(
        pathway_DB_filename,
        names=["pathway id", "num genes", "genes"],
        sep="\t",
        header=None,
    )

    # Drop extra column
    pa_pathway_DB.drop(columns=["num genes"], inplace=True)

    # Make genes tab-separated
    pa_pathway_DB["genes"] = pa_pathway_DB["genes"].str.split(";").str.join("\t")

    # Need to temporarily write data to file in order
    # to remove extra '\'
    tmp_filename = os.path.join(local_dir, "pa_pathway_DB_tmp_filename.gmt")
    pa_pathway_DB.to_csv(
        tmp_filename,
        quoting=csv.QUOTE_NONE,
        escapechar="\\",
        index=False,
        header=False,
        sep="\t",
    )
    with open(tmp_filename, "r") as ihf:
        tmp = ihf.read()
    with open(out_filename, "w") as ohf:
        ohf.write(tmp.replace("\\", ""))


def process_samples_for_limma(
    expression_filename, process_metadata_filename, grp_metadata_filename,
):
    """
    This function processes samples in the template and simulated
    experiments to prepare for DE analysis using DESeq.

    These processing steps includes:
    1. Removing samples that are not included in the comparison.
    These "extra" samples occur when an experiment contains multiple
    comparisons.
    2. Checks that the ordering of samples in the metadata file
    are consistent with the ordering in the gene expression data
    matrix. If the ordering is not consistent, then samples in
    the gene expression data matrix are re-ordered.

    Arguments
    ----------
    expression_filename: str
        File containing unnormalized gene expression data for 
        either template or simulated experiments

    process_metadata_filename: str
        File containing assignment for which samples to drop

    grp_metadata_filename: str
        File containing group assigments for samples to use
        for DESeq analysis

    """

    # Read data
    expression = pd.read_csv(expression_filename, sep="\t", index_col=0, header=0)
    process_metadata = pd.read_csv(
        process_metadata_filename, sep="\t", index_col=0, header=0
    )
    grp_metadata = pd.read_csv(grp_metadata_filename, sep="\t", header=0, index_col=0)

    # Get samples ids to remove
    samples_to_remove = list(
        process_metadata[process_metadata["processing"] == "drop"].index
    )

    # Remove samples
    expression = expression.drop(samples_to_remove)

    # Check ordering of sample ids is consistent between gene expression data and metadata
    metadata_sample_ids = grp_metadata.index
    expression_sample_ids = expression.index

    if metadata_sample_ids.equals(expression_sample_ids):
        print("sample ids are ordered correctly")
    else:
        # Convert gene expression ordering to be the same as
        # metadata sample ordering
        print("sample ids don't match, going to re-order gene expression samples")
        expression = expression.reindex(metadata_sample_ids)

    # Save
    expression.to_csv(expression_filename, sep="\t")


def process_samples_for_DESeq(
    expression_filename,
    process_metadata_filename,
    grp_metadata_filename,
    count_threshold=None,
    out_expression_filename=None,
):
    """
    This function processes samples in the template and simulated
    experiments to prepare for DE analysis using DESeq.

    These processing steps includes:
    1. Removing samples that are not included in the comparison.
    These "extra" samples occur when an experiment contains multiple
    comparisons.
    2. Removes genes with 0 counts across all samples
    3. (Optionally) filters genes with mean gene expression below
    some user defined threshold
    4. Case count values as integers
    5. Checks that the ordering of samples in the metadata file
    are consistent with the ordering in the gene expression data
    matrix. If the ordering is not consistent, then samples in
    the gene expression data matrix are re-ordered.

    Arguments
    ----------
    expression_filename: str
        File containing unnormalized gene expression data for 
        either template or simulated experiments

    process_metadata_filename: str
        File containing assignment for which samples to drop

    grp_metadata_filename: str
        File containing group assigments for samples to use
        for DESeq analysis

    count_threshold: int
        Remove genes that have mean count <= count_threshold
    
    """

    # Read data
    expression = pd.read_csv(expression_filename, sep="\t", index_col=0, header=0)
    process_metadata = pd.read_csv(
        process_metadata_filename, sep="\t", index_col=0, header=0
    )
    grp_metadata = pd.read_csv(grp_metadata_filename, sep="\t", header=0, index_col=0)

    # Get samples ids to remove
    samples_to_remove = list(
        process_metadata[process_metadata["processing"] == "drop"].index
    )

    # Remove samples
    expression = expression.drop(samples_to_remove)

    # Remove genes with 0 counts
    all_zero_genes = list(expression.columns[(expression == 0).all()])
    expression = expression.drop(all_zero_genes)

    assert len(list(expression.columns[(expression == 0).all()])) == 0

    # Remove genes below a certain threshold (if provided)
    if count_threshold != None:
        genes_to_keep = expression.loc[:, expression.mean() <= count_threshold].columns
        expression = expression[genes_to_keep]

    # Cast as int
    expression = expression.astype(int)

    # Check ordering of sample ids is consistent between gene expression data and metadata
    metadata_sample_ids = grp_metadata.index
    expression_sample_ids = expression.index

    if metadata_sample_ids.equals(expression_sample_ids):
        print("sample ids are ordered correctly")
    else:
        # Convert gene expression ordering to be the same as
        # metadata sample ordering
        print("sample ids don't match, going to re-order gene expression samples")
        expression = expression.reindex(metadata_sample_ids)

    # Save
    if out_expression_filename != None:
        expression.to_csv(out_expression_filename, sep="\t")
    else:
        expression.to_csv(expression_filename, sep="\t")
