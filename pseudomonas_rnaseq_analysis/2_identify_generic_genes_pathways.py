# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:generic_expression] *
#     language: python
#     name: conda-env-generic_expression-py
# ---

# # Identify generic genes and pathways
#
# This notebook is meant to identify common differentially expressed genes in the PAO1 and PA14 compendia

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

import os
import sys
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import pickle
import scipy.stats as ss
from keras.models import load_model
from rpy2.robjects import pandas2ri
from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import process, stats, ranking

pandas2ri.activate()

np.random.seed(123)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_pa14_rnaseq.tsv")
)
params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
NN_architecture = params["NN_architecture"]
num_runs = params["num_simulated"]
project_id = params["project_id"]
metadata_col_id = params["metadata_colname"]
raw_template_filename = params["raw_template_filename"]
processed_template_filename = params["processed_template_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]
scaler_filename = params["scaler_filename"]
col_to_rank_genes = params["rank_genes_by"]
logFC_name = params["DE_logFC_name"]
pvalue_name = params["DE_pvalue_name"]
latent_dim = params["latent_dim"]

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_groups.tsv"
)

# Load pickled file
scaler = pickle.load(open(scaler_filename, "rb"))

# Percentile threshold to identify generic genes
percentile_threshold = 80.0

metadata_simulate_filename = "data/metadata/SraRunTable.csv"
metadata_delimiter = ","
experiment_id_colname = "SRA_study"
# -

# Output files
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)


# ## Need to customize code from ponyo
#
# The current simulation-related function in ponyo, `get_sample_ids` assumes that the user is using one of two different metadata files (one associated with the pseudomonas compendium and another associated with recount2). The compendium dataset we are using here has a slightly different format for their metadata file.
#
# Here we are temporarily writing our own function customized for this Pa RNA-seq compendia

def get_sample_ids(
    metadata_filename, delimiter, experiment_colname, experiment_id, sample_id_colname
):
    """
    Returns sample ids (found in gene expression df) associated with
    a given list of experiment ids (found in the metadata)

    Arguments
    ----------
    metadata_filename: str
        Metadata file path. An example metadata file can be found
        here: https://github.com/greenelab/ponyo/blob/master/human_tests/data/metadata/recount2_metadata.tsv

    delimiter: str
        Delimiter for metadata file

    experiment_colname: str
        Column header that contains the experiment ids

    experiment_id: str
        Experiment id selected to retrieve sample ids for

    sample_id_colname: str
        Column header that contains sample id that maps expression data
        and metadata

    """

    # Read in metadata
    metadata = pd.read_csv(metadata_filename, header=0, sep=delimiter, index_col=None)

    # Set index column to experiment id column
    metadata.set_index(experiment_colname, inplace=True)

    # Select samples associated with experiment id
    selected_metadata = metadata.loc[experiment_id]
    sample_ids = list(selected_metadata[sample_id_colname])

    return sample_ids


# +
def shift_template_experiment_custom(
    normalized_data_filename,
    NN_architecture,
    latent_dim,
    dataset_name,
    scaler,
    metadata_filename,
    metadata_delimiter,
    experiment_id_colname,
    sample_id_colname,
    selected_experiment_id,
    local_dir,
    base_dir,
    num_runs,
):
    """
    Generate new simulated experiment using the selected_experiment_id as a template
    experiment using the same workflow as `simulate_by_latent_transform`

    This will return a file with a single simulated experiment following the workflow mentioned.
    This function can be run multiple times to generate multiple simulated experiments from a
    single selected_experiment_id.

    Arguments
    ----------
    normalized_data_filename: str
        File containing normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...

    NN_architecture: str
        Name of neural network architecture to use.
        Format 'NN_<intermediate layer>_<latent layer>'

    latent_dim: int
        The number of dimensions in the latent space

    dataset_name: str
        Name for analysis directory. Either "Human" or "Pseudomonas"

    scaler: minmax model
        Model used to transform data into a different range

    metadata_filename: str
        Metadata file path. Note: The format of this metadata file
        requires the index column to contain experiment ids.

    metadata_delimiter: str
        Delimiter for metadata file

    experiment_colname: str
        Column header that contains the experiment ids

    sample_id_colname: str
        Column header that contains sample id that maps expression data
        and metadata

    selected_experiment_id: str
        Experiment id selected as template

    local_dir: str
        Parent directory on local machine to store intermediate results

    base_dir: str
        Root directory containing analysis subdirectories

    num_runs: int
        Number of experiments to simulate

    Returns
    --------
    simulated_data_filename: str
        File containing simulated gene expression data

    """

    # Files
    NN_dir = os.path.join(base_dir, dataset_name, "models", NN_architecture)

    model_encoder_filename = glob.glob(os.path.join(NN_dir, "*_encoder_model.h5"))[0]

    weights_encoder_filename = glob.glob(os.path.join(NN_dir, "*_encoder_weights.h5"))[
        0
    ]

    model_decoder_filename = glob.glob(os.path.join(NN_dir, "*_decoder_model.h5"))[0]

    weights_decoder_filename = glob.glob(os.path.join(NN_dir, "*_decoder_weights.h5"))[
        0
    ]

    # Load saved models
    loaded_model = load_model(model_encoder_filename, compile=False)
    loaded_decode_model = load_model(model_decoder_filename, compile=False)

    loaded_model.load_weights(weights_encoder_filename)
    loaded_decode_model.load_weights(weights_decoder_filename)

    # Read data
    normalized_data = pd.read_csv(
        normalized_data_filename, header=0, sep="\t", index_col=0
    )

    # Get corresponding sample ids
    sample_ids = get_sample_ids(
        metadata_filename,
        metadata_delimiter,
        experiment_id_colname,
        selected_experiment_id,
        sample_id_colname,
    )

    # Gene expression data for selected samples
    selected_data_df = normalized_data.loc[sample_ids]

    for run in range(num_runs):
        simulated_data_decoded_df, simulated_data_encoded_df = run_shift_template(
            loaded_model,
            loaded_decode_model,
            normalized_data,
            selected_data_df,
            latent_dim,
        )

        # Un-normalize the data in order to run DE analysis downstream
        simulated_data_scaled = scaler.inverse_transform(simulated_data_decoded_df)

        simulated_data_scaled_df = pd.DataFrame(
            simulated_data_scaled,
            columns=simulated_data_decoded_df.columns,
            index=simulated_data_decoded_df.index,
        )

        # Save
        out_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            "selected_simulated_data_"
            + selected_experiment_id
            + "_"
            + str(run)
            + ".txt",
        )
        simulated_data_scaled_df.to_csv(out_filename, float_format="%.3f", sep="\t")

        out_encoded_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_encoded_data_{selected_experiment_id}_{run}.txt",
        )

        simulated_data_encoded_df.to_csv(
            out_encoded_filename, float_format="%.3f", sep="\t"
        )

    # Save template data for visualization validation
    test_filename = os.path.join(
        local_dir,
        "pseudo_experiment",
        "template_normalized_data_" + selected_experiment_id + "_test.txt",
    )
    selected_data_df.to_csv(test_filename, float_format="%.3f", sep="\t")


def run_shift_template(encoder, decoder, normalized_data, selected_data_df, latent_dim):
    """
    This function does the template shifting used in `shift_template_experiment`.

    Arguments
    ---------
    encoder: keras.models.Model
        The encoder half of the VAE. `encoder` takes in a (samples x genes) dataframe of
        gene expression data and encodes it into a latent space

    decoder: keras.models.Model
        The decoder half of the VAE. `decoder` takes a dataframe of means and standard deviations
        and uses them to simulate gene expression data close to the distribution of normalized_data

    normalized_data: pd.DataFrame
        The data to be used to train the VAE

    selected_data_df: pd.DataFrame
        The samples to be shifted in the latent space

    latent_dim: int
        The dimension of the latent space the samples will be shifted in

    Returns
    -------
    simulated_data_decoded_df: pd.DataFrame
        The simulated data created by shifting the samples in the latent space

    simulated_data_encoded_df: pd.DataFrame
        The latent means and standard deviations in the latent space used to simulate the data
    """
    # Encode selected experiment into latent space
    data_encoded = encoder.predict_on_batch(selected_data_df)
    data_encoded_df = pd.DataFrame(data_encoded, index=selected_data_df.index)

    # Get centroid of original data
    centroid = data_encoded_df.mean(axis=0)

    # Add individual vectors(centroid, sample point) to new_centroid

    # Encode original gene expression data into latent space
    data_encoded_all = encoder.predict_on_batch(normalized_data)
    data_encoded_all_df = pd.DataFrame(data_encoded_all, index=normalized_data.index)

    data_encoded_all_df.head()

    # Find a new location in the latent space by sampling from the latent space
    encoded_means = data_encoded_all_df.mean(axis=0)
    encoded_stds = data_encoded_all_df.std(axis=0)

    latent_dim = int(latent_dim)
    new_centroid = np.zeros(latent_dim)

    for j in range(latent_dim):
        new_centroid[j] = np.random.normal(encoded_means[j], encoded_stds[j])

    shift_vec_df = new_centroid - centroid

    simulated_data_encoded_df = data_encoded_df.apply(
        lambda x: x + shift_vec_df, axis=1
    )

    # Decode simulated data into raw gene space
    simulated_data_decoded = decoder.predict_on_batch(simulated_data_encoded_df)

    simulated_data_decoded_df = pd.DataFrame(
        simulated_data_decoded,
        index=simulated_data_encoded_df.index,
        columns=selected_data_df.columns,
    )

    return simulated_data_decoded_df, simulated_data_encoded_df


# -

# ### Simulate experiments using selected template experiment
# Workflow:
#
# 1. Get the gene expression data for the selected template experiment
# 2. Encode this experiment into a latent space using the trained VAE model
# 3. Linearly shift the encoded template experiment in the latent space
# 4. Decode the samples. This results in a new experiment
# 5. Repeat steps 1-4 to get multiple simulated experiments

# Simulate multiple experiments
# This step creates the following files in "<local_dir>/pseudo_experiment/" directory:
#   - selected_simulated_data_SRP012656_<n>.txt
#   - selected_simulated_encoded_data_SRP012656_<n>.txt
#   - template_normalized_data_SRP012656_test.txt
# in which "<n>" is an integer in the range of [0, num_runs-1]
os.makedirs(os.path.join(local_dir, "pseudo_experiment"), exist_ok=True)
shift_template_experiment_custom(
    normalized_compendium_filename,
    NN_architecture,
    latent_dim,
    dataset_name,
    scaler,
    metadata_simulate_filename,
    metadata_delimiter,
    experiment_id_colname,
    metadata_col_id,
    project_id,
    local_dir,
    base_dir,
    num_runs,
)

# +
simulated_filename = os.path.join(
    local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_1.txt"
)

test = pd.read_csv(simulated_filename, sep="\t", index_col=0, header=0)
# -

test.head()

# ### Process template and simulated data
#
# * Remove samples not required for comparison.
# * Make sure ordering of samples matches metadata for proper comparison

# +
if not os.path.exists(sample_id_metadata_filename):
    sample_id_metadata_filename = None

stats.process_samples_for_DESeq(
    raw_template_filename,
    metadata_filename,
    processed_template_filename,
    sample_id_metadata_filename,
)

for i in range(num_runs):
    simulated_filename = os.path.join(
        local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_{i}.txt"
    )
    stats.process_samples_for_DESeq(
        simulated_filename,
        metadata_filename,
        None,
        sample_id_metadata_filename,
    )

# +
# Quick check
template_data = pd.read_csv(
    processed_template_filename, header=0, index_col=0, sep="\t"
)

assert template_data.shape[0] == 6
# -

template_data.head()

# ### Differential expression analysis

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# get_DE_stats_DESeq(metadata_filename,
#                    project_id,
#                    processed_template_filename,
#                    "template",
#                    local_dir,
#                    "real")

# + magic_args="-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs -o num_sign_DEGs_simulated" language="R"
#
# source(paste0(base_dir,'/generic_expression_patterns_modules/DE_analysis.R'))
#
# num_sign_DEGs_simulated <- c()
#
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(
#         local_dir,
#         "pseudo_experiment/selected_simulated_data_",
#         project_id,
#         "_",
#         i,
#         ".txt",
#         sep=""
#     )
#
#     run_output <- get_DE_stats_DESeq(
#         metadata_filename,
#         project_id,
#         simulated_data_filename,
#         "simulated",
#         local_dir,
#         i
#     )
#     num_sign_DEGs_simulated <- c(num_sign_DEGs_simulated, run_output)
# }
# -

# ### Rank genes

analysis_type = "DE"
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)
template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
    template_DE_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_genes,
    logFC_name,
    pvalue_name,
)

# ### Gene summary table

# +
summary_gene_ranks = ranking.generate_summary_table(
    template_DE_stats_filename,
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir,
    "gene",
    params,
)

summary_gene_ranks.sort_values(by="Z score", ascending=False).head()
# -

# Add gene name as column to summary dataframe
summary_gene_ranks = ranking.add_pseudomonas_gene_name_col(summary_gene_ranks, base_dir)
summary_gene_ranks.sort_values(by="Z score", ascending=False).head()

summary_gene_ranks.sort_values(by="Percentile (simulated)", ascending=False).head()

# Check if there is an NaN values, there should not be
summary_gene_ranks.isna().any()

# Create `gene_summary_filename`
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")

# ## Compare gene ranking
#
# We can only compare the ranking between the PAO1 RNA-seq compendium vs GAPE, where we still see good concordance as expected.
#
# When we look for common genes, we do this based on percentiles generated by SOPHIE for both the PAO1 and PA14 compendia to be consistent.

if "pao1" in config_filename:
    # Get generic genes identified by Crow et. al.
    GAPE_filename = params["reference_gene_filename"]
    ref_gene_col = params["reference_gene_name_col"]
    ref_rank_col = params["reference_rank_col"]

    figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

    corr, shared_ranking = ranking.compare_gene_ranking(
        summary_gene_ranks, GAPE_filename, ref_gene_col, ref_rank_col, figure_filename
    )
