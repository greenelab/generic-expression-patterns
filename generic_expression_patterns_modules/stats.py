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
