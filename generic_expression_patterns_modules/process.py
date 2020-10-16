"""
Author: Alexandra Lee
Date Created: 16 June 2020

This script provide supporting functions to run analysis notebooks.
"""

import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from glob import glob
from sklearn.preprocessing import MinMaxScaler
from generic_expression_patterns_modules import calc
from ponyo import simulate_expression_data

# Data processing functions including:
# * functions to map ensembl gene ids to hgnc symbols
# * functions to remove subsets of samples
# * functions to transform data into integer for downstream DE analysis
# * functions to normalize data


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

    NOTE:
    -----
    This function is deprecated due to large memory usage: when `expression_df`
    is a large dataframe, manipulating it inside the momory becomes very slow
    (and sometimes even impossible) due to large memory consumption.

    The same functionality has been refactored into `get_renamed_columns()` and
    `map_recount2_data()` functions in this module.

    THIS FUNCTION IS KEPT AS A REFERENCE ONLY.
    """

    # Some columns are duplicates, for example:
    #   (ENSG00000223773.7,  ENSG00000223773) --> CD99P1
    #   (ENSG00000124334.17, ENSG00000124334) --> IL9R
    # We keep the first occurence of duplicated ensembl ids
    updated_mapping = gene_id_mapping.loc[
        ~gene_id_mapping.index.duplicated(keep="first")
    ]

    # Same ensembl ids are mapped to different gene symbol twice (CCL3L1, CCL3L3)
    # ENSG00000187510.7  ENSG00000187510    C12orf74
    # ENSG00000187510.7  ENSG00000187510    PLEKHG7
    # Manually mapping them based on what is found on ensembl site
    manual_mapping = {
        "ENSG00000187510.7": "PLEKHG7",
        "ENSG00000230417.11": "LINC00595",
        "ENSG00000255374.3": "TAS2R45",
        "ENSG00000276085.1": "CCL3L1",
    }

    # Apply manual mappings to `updated_mapping`
    for ensembl_id, gene_symbol in manual_mapping.items():
        updated_mapping.loc[ensembl_id].hgnc_symbol = gene_symbol

    # Remove paralogs.
    # Some ensembl ids are paralogs (for example, "geneA" and "geneA_PAR_Y").
    # They map to the same hgnc symbol. Homologous sequences are paralogous
    # if they were separated by a gene duplication event: if a gene in an
    # organism is duplicated to occupy two different positions in the same
    # genome, then the two copies are paralogous.
    updated_expression_df = expression_df.iloc[
        :, ~expression_df.columns.str.contains("PAR_Y")
    ]

    # Replace ensembl ids with gene symbol
    updated_expression_df.columns = updated_expression_df.columns.map(
        updated_mapping["hgnc_symbol"]
    )

    # Remove columns whose mapped ensembl id is an empty string
    updated_expression_df = updated_expression_df.iloc[
        :, updated_expression_df.columns != ""
    ]

    # Remove columns whose mapped ensembl id is `NaN`
    updated_expression_df = updated_expression_df.iloc[
        :, updated_expression_df.columns.notnull()
    ]

    return updated_expression_df


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


def recast_int(num_runs, local_dir, project_id):
    """
    Re-casts simulated experiment data to integer to use DESeq.

    Arguments
    ---------
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

        # Cast as int
        simulated_data = simulated_data.astype(int)

        # Save
        simulated_data.to_csv(simulated_data_file, float_format="%.5f", sep="\t")


def create_recount2_compendium(download_dir, output_filename):
    """
    Concatenate `t_data_counts.tsv` in each project directory and create the
    single recount2 commpendium file in TSV format.
    The first row in each `t_data_counts.tsv` is a header line that includes
    column names, so only the header in the first `t_data_counts.tsv` is copied
    to the output file.

    Arguments
    ---------
    download_dir: str
        dirname that hosts all downloaded projects data
    output_filename: str
        filename of output single compendium data
    """

    data_counts_filenames = glob(f"{download_dir}/*/t_data_counts.tsv")
    data_counts_filenames.sort()

    compendium_header = None
    with open(output_filename, "w") as ofh:
        for filename in data_counts_filenames:
            with open(filename) as ifh:
                file_header = ifh.readline()
                if compendium_header is None:
                    compendium_header = file_header
                    ofh.write(compendium_header)
                elif file_header != compendium_header:
                    raise Exception(f"Inconsistent header in {filename}")

                file_content = ifh.read()
                ofh.write(file_content)


def get_published_generic_genes(filename):
    """
    Get generic genes based on input filename, which could be a URL.

    Arguments
    ---------
    filename: str
        name of the file that includes published generic genes
    """

    df = pd.read_csv(filename, header=0, sep="\t")
    published_generic_genes = list(df["Gene_Name"])
    return published_generic_genes


def get_merged_gene_id_mapping(gene_id_filename, raw_ensembl_genes):
    """
    Merge genes in input gene_id file with the raw ensembl gene IDs.

    Arguments
    ---------
    gene_id_filename: str
        filename of input gene IDs;
    raw_ensembl_genes: list
        list of strings (ensembl gene IDs)

    Returns
    -------
    Mapping between ensembl ids and hgnc symbols
    """

    original_gene_id_mapping = pd.read_csv(
        gene_id_filename, header=0, sep="\t", index_col=0
    )

    # Get mapping between ensembl ids with and without version numbers.
    # The genes in `ensembl_genes` has version numbers at the end.
    ensembl_gene_ids = pd.DataFrame(
        data={
            "ensembl_version": raw_ensembl_genes,
            "ensembl_parsed": [gene_id.split(".")[0] for gene_id in raw_ensembl_genes],
        }
    )

    # Map ensembl gene ids with version number to gene_id_mapping
    merged_gene_id_mapping = pd.merge(
        original_gene_id_mapping,
        ensembl_gene_ids,
        left_on="ensembl_gene_id",
        right_on="ensembl_parsed",
        how="outer",
    )

    # Set `ensembl_version` column as the index
    merged_gene_id_mapping.set_index("ensembl_version", inplace=True)

    return merged_gene_id_mapping


def get_renamed_columns(
    raw_ensembl_ids,
    merged_gene_id_mapping,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
):
    """
    Find the new column names and corresponding column indexes.

    Arguments
    ---------
    raw_ensembl_ids:
        list of strings (ensembl gene IDs), which are columns names in
        raw recount2 data file;
    merged_gene_id_mapping: DataFrame
        merged gene ID mapping;
    manual_mapping: dict
        dict of manual mapping (key: ensembl_id, value: gene symbol)
    DE_prior_filename: str
        input filename that includes symbols of published generic genes
    shared_genes_filename: str
        name of output pickled file (a list of shared gene symbols)

    Returns
    -------
    A tuple that includes two entries. The first entry is a list
    of hgnc gene symbols (which will be the new column names in remapped
    recount2 data file; The second entry is a dict whose keys are hgnc gene
    symbols and values are lists of the corresponding indexes of columns in
    the raw recount2 data file (most lists include only one column index.)

    """

    updated_mapping = merged_gene_id_mapping.loc[
        ~merged_gene_id_mapping.index.duplicated(keep="first")
    ]
    for ensembl_id, gene_symbol in manual_mapping.items():
        updated_mapping.loc[ensembl_id].hgnc_symbol = gene_symbol

    # Build a dict that maps hgnc symbols to column indexes in raw recount2 data
    hgnc_to_cols = dict()
    for col_idx, ensembl_id in enumerate(raw_ensembl_ids):
        # Skip paralogs (whose ensembl_id includes "PAR_Y")
        if "PAR_Y" in ensembl_id:
            continue

        hgnc_symbol = updated_mapping.loc[ensembl_id].hgnc_symbol

        # Skip hgnc gene symbols that are `float` type (NaN) or empty strings
        if type(hgnc_symbol) == float or len(hgnc_symbol) == 0:
            continue

        if hgnc_symbol in hgnc_to_cols:
            hgnc_to_cols[hgnc_symbol].append(col_idx)
        else:
            hgnc_to_cols[hgnc_symbol] = [col_idx]

    our_gene_ids_hgnc = list(hgnc_to_cols.keys())

    published_generic_genes = get_published_generic_genes(DE_prior_filename)
    shared_genes_hgnc = list(
        set(our_gene_ids_hgnc).intersection(published_generic_genes)
    )

    # In Python, the order of elements in a list that is converted from a set
    # is non-deterministic, so it is sorted here to have reproducible result.
    shared_genes_hgnc.sort()

    # Pickle `shared_genes_hgnc` and save as `shared_genes_filename`
    if not os.path.exists(shared_genes_filename):
        with open(shared_genes_filename, "wb") as pkl_fh:
            pickle.dump(shared_genes_hgnc, pkl_fh)

    return (shared_genes_hgnc, hgnc_to_cols)


def map_recount2_data(
    raw_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
    new_filename,
):
    """
    Map the ensembl gene IDs in `raw_filename` to hgnc gene symbols based
    on the header line in `template_filename`, and save the new header
    and corresponding data columns to `new_filename`.
    """

    # Read the header line of `raw_filename` to get its column names:
    raw_header_df = pd.read_csv(raw_filename, header=0, sep="\t", nrows=1)
    raw_ensembl_ids = list(raw_header_df.columns)
    if raw_ensembl_ids[0] == "Unnamed: 0":
        del raw_ensembl_ids[0]

    merged_gene_id_mapping = get_merged_gene_id_mapping(
        gene_id_filename, raw_ensembl_ids
    )

    shared_genes_hgnc, hgnc_to_cols = get_renamed_columns(
        raw_ensembl_ids,
        merged_gene_id_mapping,
        manual_mapping,
        DE_prior_filename,
        shared_genes_filename,
    )

    col_indexes = list()
    for hgnc in shared_genes_hgnc:
        col_indexes += hgnc_to_cols[hgnc]

    output_cols = [""]
    for hgnc in shared_genes_hgnc:
        output_cols += [hgnc] * len(hgnc_to_cols[hgnc])
    output_header = "\t".join(output_cols) + "\n"

    with open(new_filename, "w") as ofh:
        ofh.write(output_header)
        with open(raw_filename) as ifh:
            for line_num, line in enumerate(ifh):
                if line_num == 0:
                    continue
                tokens = line.strip("\n").split("\t")
                sample_id = tokens[0].strip('"')
                input_values = tokens[1:]
                output_values = list()
                for idx in col_indexes:
                    output_values.append(input_values[idx])
                ofh.write(sample_id + "\t" + "\t".join(output_values) + "\n")


def process_raw_template_recount2(
    raw_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
    mapped_filename,
    sample_id_metadata_filename,
    processed_filename,
):
    """
    Create mapped recount2 template data file based on input raw template
    data file (`raw_filename`), drop sample rows if needed, and save updated
    template data on disk.
    """

    # Write the intermediate mapped recount2 template data file on disk
    map_recount2_data(
        raw_filename,
        gene_id_filename,
        manual_mapping,
        DE_prior_filename,
        shared_genes_filename,
        mapped_filename,
    )

    sample_ids_to_drop = set()
    if os.path.exists(sample_id_metadata_filename):
        # Read in metadata and get samples to be dropped:
        metadata = pd.read_csv(
            sample_id_metadata_filename, sep="\t", header=0, index_col=0
        )
        sample_ids_to_drop = set(metadata[metadata["processing"] == "drop"].index)

    # Write the processed recount2 template output file on disk
    with open(mapped_filename) as ifh, open(processed_filename, "w") as ofh:
        for idx, line in enumerate(ifh):
            sample_id = line.split("\t")[0]
            if idx == 0 or sample_id not in sample_ids_to_drop:
                ofh.write(line)


def process_raw_template_pseudomonas(
    processed_compendium_filename,
    project_id,
    dataset_name,
    metadata_colname,
    sample_id_metadata_filename,
    raw_template_filename,
    processed_template_filename,
):
    """
    Create processed pseudomonas template data file based on
    processed compendium file (`compendium_filename`),
    drop sample rows if needed, and save updated
    template data on disk.
    """

    # Get sample ids associated with selected project id
    sample_ids = simulate_expression_data.get_sample_ids(
        project_id, dataset_name, metadata_colname
    )

    # Get samples from experiment id
    processed_compendium = pd.read_csv(
        processed_compendium_filename, header=0, index_col=0, sep="\t"
    )
    template_data = processed_compendium.loc[sample_ids]

    template_data.to_csv(raw_template_filename, sep="\t")

    sample_ids_to_drop = set()
    if os.path.exists(sample_id_metadata_filename):
        # Read in metadata and get samples to be dropped:
        metadata = pd.read_csv(
            sample_id_metadata_filename, sep="\t", header=0, index_col=0
        )
        sample_ids_to_drop = set(metadata[metadata["processing"] == "drop"].index)

    # Write the processed pseudomonas template output file on disk
    with open(raw_template_filename) as ifh, open(
        processed_template_filename, "w"
    ) as ofh:
        for idx, line in enumerate(ifh):
            sample_id = line.split("\t")[0]
            if idx == 0 or sample_id not in sample_ids_to_drop:
                ofh.write(line)


def normalize_compendium(
    mapped_filename, normalized_filename, scaler_filename,
):
    """
    Read the mapped compendium file into memory, normalize it, and save
    both normalized compendium data and pickled scaler on disk.
    """

    # Read mapped compendium file: ~4 minutes (17 GB of RAM)
    mapped_compendium_df = pd.read_table(
        mapped_filename, header=0, sep="\t", index_col=0
    )

    # 0-1 normalize per gene
    scaler = MinMaxScaler()

    # Fitting (2 minutes, ~8 GB of RAM)
    normalized_compendium = scaler.fit_transform(mapped_compendium_df)
    normalized_compendium_df = pd.DataFrame(
        normalized_compendium,
        columns=mapped_compendium_df.columns,
        index=mapped_compendium_df.index,
    )

    # Save normalized data on disk: ~17.5 minutes
    normalized_compendium_df.to_csv(normalized_filename, float_format="%.3f", sep="\t")

    # Pickle `scaler` as `scaler_filename` on disk
    with open(scaler_filename, "wb") as pkl_fh:
        pickle.dump(scaler, pkl_fh)


def process_raw_compendium_pseudomonas(
    raw_filename, processed_filename, normalized_filename, scaler_filename,
):
    """
    Create processed pseudomonas compendium data file based on raw compendium
    data file (`raw_filename`), and normalize the processed compendium.
    """

    # Create processed pseudomonas compendium data file
    raw_compendium = pd.read_csv(raw_filename, header=0, index_col=0, sep="\t")

    if raw_compendium.shape != (950, 5549):
        processed_compendium = raw_compendium.T
    else:
        processed_compendium = raw_compendium

    assert processed_compendium.shape == (950, 5549)

    # Save transformed compendium data
    processed_compendium.to_csv(processed_filename, sep="\t")

    # Normalize processed pseudomonas compendium data
    normalize_compendium(processed_filename, normalized_filename, scaler_filename)


# Functions to format intermediate data files to prepare to compare gene/pathway
# ranking:
# * functions to concatenate simulated data results
# * functions to get absolute value of test statistics to use for ranking
# * functions to generate summary data files
# * functions to scale ranking


def concat_simulated_data(local_dir, num_runs, project_id, data_type):
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
    data_type: str
        Either 'DE' or 'GSEA'

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
    This function takes the absolute value of columns=[`logFC`, `t`].
    For ranking genes, we only care about the magnitude of the change for
    the logFC and t statistic, but not the direction.

    The ranking for each gene will be based on the mean absolute value of either
    logFC or t statistic, depending on the user selection
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
    if "adj.P.Val" in template_simulated_DE_stats.columns:
        median_pval_simulated = template_simulated_DE_stats[("adj.P.Val", "median")]
        col_name = "adj.P.Val"
    else:
        median_pval_simulated = template_simulated_DE_stats[("padj", "median")]
        col_name = "padj"
    mean_test_simulated = template_simulated_DE_stats[(col_to_rank, "mean")]
    std_test_simulated = template_simulated_DE_stats[(col_to_rank, "std")]
    count_simulated = template_simulated_DE_stats[(col_to_rank, "count")]
    rank_simulated = template_simulated_DE_stats[("ranking", "")]

    summary = pd.DataFrame(
        data={
            "Gene ID": template_simulated_DE_stats.index,
            "Adj P-value (Real)": template_simulated_DE_stats[col_name],
            "Rank (Real)": template_simulated_DE_stats["ranking"],
            "Test statistic (Real)": template_simulated_DE_stats[col_to_rank],
            "Median adj p-value (simulated)": median_pval_simulated,
            "Rank (simulated)": rank_simulated,
            "Mean test statistic (simulated)": mean_test_simulated,
            "Std deviation (simulated)": std_test_simulated,
            "Number of experiments (simulated)": count_simulated,
        }
    )
    summary["abs(Z score)"] = (
        abs(
            summary["Test statistic (Real)"]
            - summary["Mean test statistic (simulated)"]
        )
    ) / summary["Std deviation (simulated)"]

    # Save file
    summary_file = os.path.join(local_dir, "gene_summary_table_" + col_to_rank + ".tsv")

    summary.to_csv(summary_file, float_format="%.5f", sep="\t")

    return summary


def merge_ranks_to_compare(
    your_summary_ranks_df, reference_ranks_file, reference_name_col, reference_rank_col,
):
    """
    Given dataframes of your ranking of genes or pathways
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


def get_shared_rank_scaled(
    summary_df, reference_filename, ref_gene_col, ref_rank_col, data_type
):
    """
    Returns shared rank scaled dataframe and correlation values based on
    input `summary_df` and other parameters.

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
    data_type: str
        Either 'DE' or 'GSEA'

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

    # Note: These lowly expressed genes were not pre-filtered before DESeq
    # (Micheal Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of
    # dispersion is difficult for genes with very, very low average counts. See the methods.
    # However, it doesn't really effect the outcome because these genes have almost no power for
    # detecting differential expression. Effects runtime though.
    shared_rank_scaled_df = shared_rank_scaled_df[
        ~shared_rank_scaled_df["Rank (simulated)"].isna()
    ]

    # Get correlation
    r, p, ci_low, ci_high = calc.spearman_ci(
        0.95, shared_rank_scaled_df, 1000, data_type
    )

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
        "Our preliminary method", "DE prior (Crow et. al. 2019)", fontsize=14
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
        color="slateblue",
        s=100,
    )

    fig.set_xlabel("Our preliminary method", fontsize=14)
    fig.set_ylabel("Powers et. al. 2018", fontsize=14)

    fig.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )

    return correlations


def concat_simulated_data_columns(local_dir, num_runs, project_id, data_type):
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
    data_type: str
        Either 'DE' or 'GSEA'

    Returns
    -------
    Dataframe containing all simulated experiments concatenated together

    """

    # Only "DE" and "GSEA" data types are supported.
    if data_type not in ["DE", "GSEA"]:
        raise Exception(f"Unknown data_type: {data_type}")

    simulated_stats_all = pd.DataFrame()
    for i in range(num_runs):
        simulated_stats_filename = os.path.join(
            local_dir,
            f"{data_type}_stats",
            f"{data_type}_stats_simulated_data_{project_id}_{i}.txt",
        )
        # Read results
        simulated_stats = pd.read_csv(
            simulated_stats_filename, header=0, sep="\t", index_col=0
        )

        # Concatenate df
        simulated_stats_all = pd.concat(
            [simulated_stats_all, simulated_stats["NES"]], axis=1
        )

    simulated_stats_all.index = simulated_stats.index
    return simulated_stats_all


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


# Functions related to visualizing trends in generic
# genes/pathways found
# * functions to generate summary dataframes
# * functions to plot trends
# * functions to compare groups of genes


def merge_abs_raw_dfs(abs_df, raw_df, condition):
    """
    This function merges and returns dataframe containing
    summary gene results using absolute value of the test 
    statistic and raw test statistic values.

    Arguments
    ---------
    abs_df: df
        Summary df using absolute value of test statistic
    raw_df: df
        Summary df using raw value of test statistic
    condition: str
        Condition from E-GEOD-33245. Either '1v2', '1v3', '1v4' or '1v5'
    """
    merged_df = abs_df.merge(
        raw_df,
        left_on="Gene ID",
        right_on="Gene ID",
        suffixes=[f"_grp_{condition}", f"_grp_{condition}_raw"],
    )

    return merged_df


def merge_two_conditions_df(
    merged_condition_1_df, merged_condition_2_df, condition_1, condition_2
):
    """
    This function merges and returns summary dataframes across two conditions to
    compare trends. For example, merge summary dataframes between 1v2 and 1v3.

    Arguments
    ---------
    merged_condition_1_df: df
        df of results for one of the E-GEOD-33245 conditions ('1v2', '1v3', '1v4' or '1v5')
        returned from `merge_abs_raw_dfs`
    merged_condition_2_df: df
        df of results for another one of the E-GEOD-33245 conditions ('1v2', '1v3', '1v4' or '1v5')
        returned from `merge_abs_raw_dfs`
    condition_1: str
        Condition from E-GEOD-33245 associated with 'merged_condition_1_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    condition_2: str
        Condition from E-GEOD-33245 associated with 'merged_condition_2_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    """
    merged_all_df = merged_condition_1_df.merge(
        merged_condition_2_df, left_on="Gene ID", right_on="Gene ID"
    )
    merged_all_df["max Z score"] = (
        merged_all_df[
            [f"abs(Z score)_grp_{condition_1}", f"abs(Z score)_grp_{condition_2}"]
        ]
        .abs()
        .max(axis=1)
    )
    merged_all_df["Gene ID Name"] = (
        merged_all_df["Gene ID"]
        + " "
        + merged_all_df[f"Gene Name_grp_{condition_1}"].fillna("")
    )

    merged_df = merged_all_df[
        [
            "Gene ID",
            "Gene ID Name",
            f"Test statistic (Real)_grp_{condition_1}",
            f"Test statistic (Real)_grp_{condition_1}_raw",
            f"Adj P-value (Real)_grp_{condition_1}",
            f"Mean test statistic (simulated)_grp_{condition_1}",
            f"Std deviation (simulated)_grp_{condition_1}",
            f"Median adj p-value (simulated)_grp_{condition_1}",
            f"Test statistic (Real)_grp_{condition_2}",
            f"Test statistic (Real)_grp_{condition_2}_raw",
            f"Adj P-value (Real)_grp_{condition_2}",
            f"Mean test statistic (simulated)_grp_{condition_2}",
            f"Std deviation (simulated)_grp_{condition_2}",
            f"Median adj p-value (simulated)_grp_{condition_2}",
            f"abs(Z score)_grp_{condition_1}",
            f"abs(Z score)_grp_{condition_2}",
            "max Z score",
        ]
    ]
    return merged_df


def plot_two_conditions(merged_df, condition_1, condition_2, xlabel, ylabel):
    """
    This function plots scatterplot comparing trends across two
    conditions

    Arguments
    ---------
    merged_df: df
        Merged df containing results for two conditions of E-GEOD-33245.
        Created from `merge_two_conditions_df`
    condition_1:condition_1: str
        Condition from E-GEOD-33245 associated with 'merged_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    condition_2: str
        Condition from E-GEOD-33245 associated with 'merged_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    xlabel: str
        Label to describe condition_1
    ylabel: str
        Label to describe condition_2

    """
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))
    cmap = sns.cubehelix_palette(start=2.8, rot=0.1, as_cmap=True)

    fig_abs = sns.scatterplot(
        data=merged_df,
        x=f"Test statistic (Real)_grp_{condition_1}",
        y=f"Test statistic (Real)_grp_{condition_2}",
        hue="max Z score",
        size="max Z score",
        linewidth=0,
        alpha=0.7,
        ax=axes[0],
        palette=cmap,
    )
    fig_abs.plot([0, 4], [0, 4], "--k")

    fig_raw = sns.scatterplot(
        data=merged_df,
        x=f"Test statistic (Real)_grp_{condition_1}_raw",
        y=f"Test statistic (Real)_grp_{condition_2}_raw",
        hue="max Z score",
        size="max Z score",
        linewidth=0,
        alpha=0.7,
        ax=axes[1],
        palette=cmap,
    )
    fig_raw.plot([-4, 4], [-4, 4], "--k")

    # Add labels
    fig.suptitle(f"({xlabel}) vs ({ylabel})", fontsize=16)
    fig.text(0.5, 0.04, xlabel, ha="center", va="center")
    fig.text(0.06, 0.5, ylabel, ha="center", va="center", rotation="vertical")
    axes[0].set_title("using abs(log$_2$ Fold Change)")
    axes[1].set_title("using log$_2$ Fold Change")
    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[0].set_ylabel("")
    axes[1].set_ylabel("")
    print(fig)


def get_and_save_DEG_lists(
    merged_one_condition_df, condition, p_threshold, z_threshold
):
    """
    Get list of DEGs using traditional criteria (log2FC and p-value)
    and using z-score cutoff. Return different combinations of gene
    lists.

    Arguments
    ---------
    merged_one_condition_df: df
        df of results for one of the E-GEOD-33245 conditions ('1v2', '1v3', '1v4' or '1v5')
        returned from `merge_abs_raw_dfs`
    condition: str
        Condition from E-GEOD-33245 associated with 'merged_one_condition_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    """
    # Get DEGs using traditional criteria
    degs_traditional = list(
        (
            merged_one_condition_df[
                (merged_one_condition_df[f"Test statistic (Real)_grp_{condition}"] > 1)
                & (
                    merged_one_condition_df[f"Adj P-value (Real)_grp_{condition}"]
                    < p_threshold
                )
            ]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of DEGs using traditional criteria: {len(degs_traditional)}")

    # Get predicted specific DEGs using z-score cutoff
    degs_specific = list(
        (
            merged_one_condition_df[
                (merged_one_condition_df[f"Test statistic (Real)_grp_{condition}"] > 1)
                & (
                    merged_one_condition_df[f"abs(Z score)_grp_{condition}"].abs()
                    > z_threshold
                )
            ]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of specific DEGs using z-score: {len(degs_specific)}")

    # Get predicted generic DEGs using z-score cutoff
    # Z-score cutoff was found by calculating the score
    # whose invnorm(0.05/5549). Here we are using a p-value = 0.05
    # with a Bonferroni correction for 5549 tests, which are
    # the number of P. aeruginosa genes
    degs_generic = list(
        (
            merged_one_condition_df[
                (merged_one_condition_df[f"Test statistic (Real)_grp_{condition}"] > 1)
                & (
                    merged_one_condition_df[f"abs(Z score)_grp_{condition}"].abs()
                    < z_threshold
                )
            ]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of generic DEGs using z-score: {len(degs_generic)}")

    # Get intersection of DEGs using traditional and z-score criteria
    degs_intersect = list(set(degs_traditional).intersection(degs_specific))
    print(
        f"No. of traditional DEGs that are specific by z-score criteria: {len(degs_intersect)}"
    )

    # Get specific DEGs that were NOT found using traditional criteria
    degs_diff = list(set(degs_specific).difference(degs_intersect))
    print(
        f"No. of specific DEGs that were not found by traditional criteria: {len(degs_diff)}"
    )

    # Get intersection of DEGs using traditional and z-score criteria
    degs_intersect_generic = list(set(degs_traditional).intersection(degs_generic))
    print(
        f"No. of traditional DEGs that are generic by z-score criteria: {len(degs_intersect_generic)}"
    )

    # Save list of genes that interesect and those that do not
    merged_one_condition_df["Gene ID Name"] = (
        merged_one_condition_df["Gene ID"]
        + " "
        + merged_one_condition_df[f"Gene Name_grp_{condition}"].fillna("")
    )

    # Set `Gene ID` as index
    merged_one_condition_df.set_index("Gene ID", inplace=True)

    gene_id_names_intersect = merged_one_condition_df.loc[
        degs_intersect, "Gene ID Name"
    ]
    gene_id_names_diff = merged_one_condition_df.loc[degs_diff, "Gene ID Name"]
    gene_id_names_generic = merged_one_condition_df.loc[degs_generic, "Gene ID Name"]

    gene_lists_df = pd.DataFrame(
        {
            "Traditional + specific DEGs": gene_id_names_intersect,
            "Specific only DEGs": gene_id_names_diff,
            "Generic DEGs": gene_id_names_generic,
        }
    )

    return (
        gene_lists_df,
        degs_traditional,
        degs_specific,
        degs_generic,
        degs_intersect,
        degs_intersect_generic,
        degs_diff,
    )


def plot_volcanos(
    degs_intersect, degs_diff, merged_one_condition_df, condition, fig_title
):
    """
    Make volcano plots based on one condition from E-GEOD-33245. Color genes
    by gene lists created from `get_and_save_DEG_lists`

    Arguments
    ---------
    degs_intersect: list
        List of genes that were found to be DE using traditional criteria
        and were found to have a high z-score (specificity)
    degs_diff: list
        List of genes that were found to have a high log2 fold change and
        high z-score but were not found to be DE using traditional criteria
    merged_one_condition_df: df
        df of results for one of the E-GEOD-33245 conditions ('1v2', '1v3', '1v4' or '1v5')
        returned from `merge_abs_raw_dfs`
    condition: str
        Condition from E-GEOD-33245 associated with 'merged_one_condition_df'.
        Either '1v2', '1v3', '1v4' or '1v5'
    fig_title: str
        Title to describe condition
    """
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))

    # Add columns for plotting
    merged_one_condition_df["FDR adjusted p-value plot"] = -np.log10(
        merged_one_condition_df[f"Adj P-value (Real)_grp_{condition}"]
    )
    merged_one_condition_df["gene group"] = "none"
    merged_one_condition_df.loc[
        degs_intersect, "gene group"
    ] = "traditional + specific DEGs"
    merged_one_condition_df.loc[degs_diff, "gene group"] = "only specific DEGs"

    colors = ["lightgrey", "red", "blue"]
    # Plot: log2FC vs p-value (traditional criteria)
    sns.scatterplot(
        data=merged_one_condition_df,
        x=f"Test statistic (Real)_grp_{condition}_raw",
        y="FDR adjusted p-value plot",
        hue="gene group",
        hue_order=["none", "traditional + specific DEGs", "only specific DEGs"],
        style="gene group",
        markers={
            "none": ".",
            "traditional + specific DEGs": "o",
            "only specific DEGs": "o",
        },
        palette=colors,
        linewidth=0,
        alpha=0.5,
        ax=axes[0],
    )

    # Plot: log2FC vs z-score
    sns.scatterplot(
        data=merged_one_condition_df,
        x=f"Test statistic (Real)_grp_{condition}_raw",
        y=f"abs(Z score)_grp_{condition}",
        hue="gene group",
        hue_order=["none", "traditional + specific DEGs", "only specific DEGs"],
        style="gene group",
        markers={
            "none": ".",
            "traditional + specific DEGs": "o",
            "only specific DEGs": "o",
        },
        palette=colors,
        linewidth=0,
        alpha=0.5,
        ax=axes[1],
    )

    # Plot: z-score vs p-value
    sns.scatterplot(
        data=merged_one_condition_df,
        x=f"abs(Z score)_grp_{condition}",
        y="FDR adjusted p-value plot",
        hue="gene group",
        hue_order=["none", "traditional + specific DEGs", "only specific DEGs"],
        style="gene group",
        markers={
            "none": ".",
            "traditional + specific DEGs": "o",
            "only specific DEGs": "o",
        },
        palette=colors,
        linewidth=0,
        alpha=0.5,
        ax=axes[2],
    )

    # Add labels
    fig.suptitle(fig_title, fontsize=16)
    axes[0].set_xlabel("log$_2$ Fold Change")
    axes[1].set_xlabel("log$_2$ Fold Change")
    axes[2].set_xlabel("Z-score")
    axes[0].set_ylabel("FDR adjusted p-value")
    axes[1].set_ylabel("Z-score")
    axes[2].set_ylabel("FDR adjusted p-value")
    axes[0].set_title("log$_2$ Fold Change vs p-value")
    axes[1].set_title("log$_2$ Fold Change vs z-score")
    axes[2].set_title("z-score vs p-value")
    print(fig)


def plot_venn(degs_traditional, degs_specific, degs_generic):
    """
    Create venn diagram to compare the genes that were found
    to be DE using traditional criteria vs genes that are
    specific (i.e. high z-score) or generic (i.e. low z-score)

    Arguments
    ---------
    degs_traditional: list
        List of genes found to pass traditional DE criteria
        (log2FC > 1 and FDR adjusted p-value < 0.05).
    degs_specific: list
        List of genes that were found to have log2 FC > 1
        and z-score > 4.44
    degs_generic: list
        List of genes that were found to have log2 FC > 1
        and z-score < 4.44
    """
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(15, 4))

    venn2(
        [set(degs_traditional), set(degs_specific)],
        set_labels=("Traditional", "Specific"),
        ax=axes[0],
    )

    venn2(
        [set(degs_traditional), set(degs_generic)],
        set_labels=("Traditional", "Generic"),
        ax=axes[1],
    )
