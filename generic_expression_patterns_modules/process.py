"""
Author: Alexandra Lee
Date Created: 16 June 2020

These scripts provide supporting functions to run analysis notebooks.
These scripts include: replacing ensembl gene ids with hgnc symbols.
"""

import pandas as pd
import os
import numpy as np
import sklearn
from glob import glob
from sklearn.preprocessing import MinMaxScaler


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
    summary["Z score"] = (
        summary["Test statistic (Real)"] - summary["Mean test statistic (simulated)"]
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
        the dirname that hosts all downloaded projects data
    output_filename: str
        the filename of the output single compendium data
    """

    data_counts_filenames = glob(f"{download_dir}/*/t_data_counts.tsv")
    data_counts_filenames.sort()

    compendium_header = None
    with open(output_filename, 'w') as ofh:
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
    published_generic_genes = list(df['Gene_Name'])
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
        gene_id_filename, header=0, sep='\t', index_col=0
    )

    # Get mapping between ensembl ids with and without version numbers.
    # The genes in `ensembl_genes` has version numbers at the end.
    ensembl_gene_ids = pd.DataFrame(
        data={
            'ensembl_version': raw_ensembl_genes,
            'ensembl_parsed': [gene_id.split('.')[0] for gene_id in raw_ensembl_genes]
        }
    )

    # Map ensembl gene ids with version number to gene_id_mapping
    merged_gene_id_mapping = pd.merge(
        original_gene_id_mapping,
        ensembl_gene_ids,
        left_on='ensembl_gene_id',
        right_on='ensembl_parsed',
        how='outer'
    )

    # Set `ensembl_version` column as the index
    merged_gene_id_mapping.set_index('ensembl_version', inplace=True)

    return merged_gene_id_mapping


def get_renamed_columns(
        raw_ensembl_ids,
        merged_gene_id_mapping,
        manual_mapping,
        DE_prior_filename
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

    return (shared_genes_hgnc, hgnc_to_cols)


def map_recount2_data(
        raw_filename,
        gene_id_filename,
        manual_mapping,
        DE_prior_filename,
        new_filename
):
    """
    Map the ensembl gene IDs in `raw_filename` to hgnc gene symbols based
    on the header line in `template_filename`, and save the new  header
    and corresponding data columns to `new_filename`.
    """

    # Read the header line of `raw_filename` to get its column names:
    raw_header_df = pd.read_csv(raw_filename, header=0, sep="\t", nrows=1)
    raw_ensembl_ids = list(raw_header_df.columns)
    if raw_ensembl_ids[0] == 'Unnamed: 0':
        del raw_ensembl_ids[0]

    merged_gene_id_mapping = get_merged_gene_id_mapping(
        gene_id_filename, raw_ensembl_ids
    )

    shared_genes_hgnc, hgnc_to_cols = get_renamed_columns(
        raw_ensembl_ids, merged_gene_id_mapping, manual_mapping, DE_prior_filename
    )

    col_indexes = list()
    for hgnc in shared_genes_hgnc:
        col_indexes += hgnc_to_cols[hgnc]

    output_cols = [""]
    for hgnc in shared_genes_hgnc:
        output_cols += [hgnc] * len(hgnc_to_cols[hgnc])
    output_header = "\t".join(output_cols) + "\n"

    with open(new_filename, 'w') as ofh:
        ofh.write(output_header)
        with open(raw_filename) as ifh:
            for line_num, line in enumerate(ifh):
                if line_num == 0:
                    continue
                tokens = line.strip('\n').split('\t')
                sample_id = tokens[0].strip('"')
                input_values = tokens[1:]
                output_values = list()
                for idx in col_indexes:
                    output_values.append(input_values[idx])
                ofh.write(sample_id + "\t" + "\t".join(output_values) + "\n")


def process_raw_template(
        raw_filename,
        gene_id_filename,
        manual_mapping,
        DE_prior_filename,
        mapped_filename,
        sample_id_metadata_filename,
        processed_filename
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
        mapped_filename
    )

    sample_ids_to_drop = set()
    if os.path.exists(sample_id_metadata_filename):
        # Read in metadata and get samples to be dropped:
        metadata = pd.read_csv(
            sample_id_metadata_filename, sep='\t', header=0, index_col=0
        )
        sample_ids_to_drop = set(metadata[metadata["processing"] == "drop"].index)

    # Write the processed recount2 template output file on disk
    with open(mapped_filename) as ifh, open(processed_filename, "w") as ofh:
        for idx, line in enumerate(ifh):
            sample_id = line.split('\t')[0]
            if idx == 0 or sample_id not in sample_ids_to_drop:
                ofh.write(line)


def normalize_compendium(mapped_filename, normalized_filename):
    """
    Read the mapped compendium file into memor, normalize it, then save
    the normalized compendium as a tsv file, and pickle the scaler.
    """

    # Read mapped compendium file: ~4 minutes (17 GB of RAM)
    mapped_compendium_df = pd.read_table(
        mapped_filename,
        header=0,
        sep='\t',
        index_col=0
    )

    # 0-1 normalize per gene
    scaler = MinMaxScaler()

    # Fitting (2 minutes, ~8 GB of RAM)
    normalized_compendium = scaler.fit_transform(mapped_compendium_df)
    normalized_compendium_df = pd.DataFrame(
        normalized_compendium,
        columns=mapped_compendium_df.columns,
        index=mapped_compendium_df.index
    )

    # Save normalized data on disk: ~17.5 minutes
    normalized_compendium_df.to_csv(
        normalized_filename, float_format='%.3f', sep='\t'
    )


def process_raw_compendium(
        raw_filename,
        gene_id_filename,
        manual_mapping,
        DE_prior_filename,
        mapped_filename,
        normalized_filename
):
    """
    Create mapped recount2 compendium data file based on raw compendium
    data file (`raw_filename`), then normalize the mapped compendium and
    save the data on disk.
    """

    # Create mapped recount2 compendium data file
    map_recount2_data(
        raw_filename,
        gene_id_filename,
        manual_mapping,
        DE_prior_filename,
        mapped_filename
    )

    # Normalize mapped recount2 compendium data and save it on disk
    normalize_compendium(mapped_filename, normalized_filename)
