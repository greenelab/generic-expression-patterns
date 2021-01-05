"""
Author: Alexandra Lee
Date Created: 9 December 2020

This script provides supporting functions to run analysis to identify
generic and specific genes for a new experiment of interest using an already trained VAE model.
"""

import os
import pickle
import numpy as np
import pandas as pd
import glob
from keras.models import load_model
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def transpose_save(input_gene_expression_filename, output_gene_expression_filename):
    """
    This function transposes the gene expression matrix in `input_gene_expression_filename`
    and saves the newly transposed matrix to `output_gene_expression_filename`
    """
    experiment = pd.read_csv(
        input_gene_expression_filename, sep="\t", index_col=0, header=0
    ).T

    experiment.to_csv(output_gene_expression_filename, sep="\t")


def compare_match_features(template_filename, compendium_filename):
    """
    This function checks that the feature space matches between
    template experiment and VAE model.  
    (i.e. ensure genes in template and VAE model are the same).
    
    If there are differences this function does the following:
    If a gene is present in template experiment but not in the VAE model, then drop gene
    If a gene is present in VAE model but not in the template experiment, 
    then add gene to template experiment with median gene expression value
    
    template_filename: str
        File containing template gene expression data. Expect matrix of dimension: sample x gene
        
    compendium_filename: str
        File containing un-normalized compendium gene expression data. 
        Gene ids are either using PA#### (P. aeruginosa)
        or using HGNC symbols (Human)
        
    """
    # Read template experiment
    template_experiment = pd.read_csv(
        template_filename, sep="\t", index_col=0, header=0
    )

    print(template_experiment.shape)

    # Read compendium
    compendium = pd.read_csv(compendium_filename, sep="\t", index_col=0, header=0)

    print(compendium.shape)

    # Check if genes are shared:
    template_genes = template_experiment.columns
    compendium_genes = compendium.columns

    # If a gene is present in template experiment but not in the VAE model, then drop gene
    # If a gene is present in VAE model but not in the template experiment,
    # then add gene to template experiment with median gene expression value
    only_template_genes = list(set(template_genes).difference(compendium_genes))
    only_compendium_genes = list(set(compendium_genes).difference(template_genes))

    tmp_template_experiment = template_experiment.drop(columns=only_template_genes)

    # Get median gene expression for only_compendium_genes
    # Use mapped_compendium for this to get counts
    median_gene_expression = compendium[only_compendium_genes].median().to_dict()
    tmp2_template_experiment = tmp_template_experiment.assign(**median_gene_expression)

    assert len(tmp2_template_experiment.columns) == len(compendium.columns)

    # sort template experiment columns to be in the same order as the compendium
    mapped_template_experiment = tmp2_template_experiment[compendium.columns]

    mapped_template_experiment.to_csv(template_filename, sep="\t")

    return mapped_template_experiment


def normalize_template_experiment(mapped_template_experiment, scaler_filename):
    """
    This function normalizes the template experiment to be within
    0-1 range, using the same scaler transform that was used to
    0-1 scale the training compendium.

    mapped_template_experiment: df
        Dataframe of template experiment after mapping gene ids        
    scaler_filename: str
        Filename containing picked scaler transform used to normalize compendium data
    """
    # Load pickled file
    with open(scaler_filename, "rb") as scaler_fh:
        scaler = pickle.load(scaler_fh)

    processed_template_experiment = scaler.transform(mapped_template_experiment)

    processed_template_experiment_df = pd.DataFrame(
        processed_template_experiment,
        columns=mapped_template_experiment.columns,
        index=mapped_template_experiment.index,
    )

    return processed_template_experiment_df


def process_template_experiment(
    template_filename,
    compendium_filename,
    scaler_filename,
    mapped_template_filename,
    processed_template_filename,
):
    """
    This function processes the template experiment to prepare for
    simulating new data. Specifically this function does the following:
    
    1. Compares and maps the template feature space to the compendium 
    feature space using `compare_match_features()`
    2. Normalizes the template experiment to be in the same scale
    as the compendium dataset using `normalize_template_experiment()`

    Arguments
    ----------
    template_filename: str
        File containing template gene expression data. Expect matrix of dimension: sample x gene        
    compendium_filename: str
        File containing un-normalized compendium gene expression data. 
        Gene ids are either using PA#### (P. aeruginosa)
        or using HGNC symbols (Human)
    scaler_filename: str
        Filename containing pickled scaler transform used to normalize compendium data
    mapped_filename: str
        Filename containing the template data where genes are mapped to compendium data.
    processed_filename: str
        Filename containing the template normalized data. This data can now be
        encoded into the learned latent space.
    """

    # Compare and map genes from template experiment to
    # compendium dataset
    mapped_template_experiment = compare_match_features(
        template_filename, compendium_filename
    )

    normalized_template_experiment = normalize_template_experiment(
        mapped_template_experiment, scaler_filename
    )

    # Save
    mapped_template_experiment.to_csv(mapped_template_filename, sep="\t")
    normalized_template_experiment.to_csv(processed_template_filename, sep="\t")


def embed_shift_template_experiment(
    normalized_data,
    template_experiment,
    vae_model_dir,
    selected_experiment_id,
    scaler_filename,
    local_dir,
    latent_dim,
    run,
):
    """
    Generate new simulated experiment using the selected_experiment_id as a template
    experiment and linearly shift template experiment to different locations of the
    latent space to create new experiment. This workflow is similar to `simulate_by_latent_transform`

    This will return a file with a single simulated experiment following the workflow mentioned.
    This function can be run multiple times to generate multiple simulated experiments from a
    single selected_experiment_id.

    Arguments
    ----------
    normalized_data: df
        Normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...

    scaler: minmax model
        Model used to transform data into a different range

    local_dir: str
        Parent directory on local machine to store intermediate results

    base_dir: str
        Root directory containing analysis subdirectories

    run: int
        Simulation run

    Returns
    --------
    simulated_data_filename: str
        File containing simulated gene expression data

    """

    # Files
    NN_dir = vae_model_dir

    model_encoder_filename = glob.glob(os.path.join(NN_dir, "*_encoder_model.h5"))[0]

    weights_encoder_filename = glob.glob(os.path.join(NN_dir, "*_encoder_weights.h5"))[
        0
    ]

    model_decoder_filename = glob.glob(os.path.join(NN_dir, "*_decoder_model.h5"))[0]

    weights_decoder_filename = glob.glob(os.path.join(NN_dir, "*_decoder_weights.h5"))[
        0
    ]

    # Load pickled file
    with open(scaler_filename, "rb") as scaler_fh:
        scaler = pickle.load(scaler_fh)

    # Load saved models
    loaded_model = load_model(model_encoder_filename, compile=False)
    loaded_decode_model = load_model(model_decoder_filename, compile=False)

    loaded_model.load_weights(weights_encoder_filename)
    loaded_decode_model.load_weights(weights_decoder_filename)

    # Gene expression data for template experiment
    selected_data_df = template_experiment

    # Encode selected experiment into latent space
    data_encoded = loaded_model.predict_on_batch(selected_data_df)
    data_encoded_df = pd.DataFrame(data_encoded, index=selected_data_df.index)

    # Get centroid of original data
    centroid = data_encoded_df.mean(axis=0)

    # Add individual vectors(centroid, sample point) to new_centroid

    # Encode original gene expression data into latent space
    data_encoded_all = loaded_model.predict_on_batch(normalized_data)
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
    # print(shift_vec_df)

    simulated_data_encoded_df = data_encoded_df.apply(
        lambda x: x + shift_vec_df, axis=1
    )

    # Decode simulated data into raw gene space
    simulated_data_decoded = loaded_decode_model.predict_on_batch(
        simulated_data_encoded_df
    )

    simulated_data_decoded_df = pd.DataFrame(
        simulated_data_decoded,
        index=simulated_data_encoded_df.index,
        columns=selected_data_df.columns,
    )

    # Un-normalize the data in order to run DE analysis downstream
    simulated_data_scaled = scaler.inverse_transform(simulated_data_decoded_df)

    simulated_data_scaled_df = pd.DataFrame(
        simulated_data_scaled,
        columns=simulated_data_decoded_df.columns,
        index=simulated_data_decoded_df.index,
    )

    # Save template data for visualization validation
    test_filename = os.path.join(
        local_dir,
        "pseudo_experiment",
        "template_normalized_data_" + selected_experiment_id + "_test.txt",
    )

    selected_data_df.to_csv(test_filename, float_format="%.3f", sep="\t")

    # Save
    out_filename = os.path.join(
        local_dir,
        "pseudo_experiment",
        "selected_simulated_data_" + selected_experiment_id + "_" + str(run) + ".txt",
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


def get_and_save_DEG_lists(summary_df, test_statistic, p_threshold, z_threshold):
    """
    Get list of DEGs using traditional criteria (log2FC and p-value)
    and using z-score cutoff. Return different combinations of gene
    lists.

    Arguments
    ---------
    summary_df: df
        df of results for one of the E-GEOD-33245 conditions ('1v2', '1v3', '1v4' or '1v5')
        returned from `merge_abs_raw_dfs`
    """
    # Get DEGs using traditional criteria
    degs_traditional = list(
        (
            summary_df[
                (summary_df[f"abs({test_statistic}) (Real)"] > 1)
                & (summary_df[f"Adj P-value (Real)"] < p_threshold)
            ]
            .set_index("Gene ID")
            .index
        )
    )
    print(f"No. of DEGs using traditional criteria: {len(degs_traditional)}")

    # Get predicted specific DEGs using z-score cutoff
    degs_specific = list(
        (
            summary_df[
                (summary_df[f"abs({test_statistic}) (Real)"] > 1)
                & (summary_df[f"abs(Z score)"].abs() > z_threshold)
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
            summary_df[
                (summary_df[f"abs({test_statistic}) (Real)"] > 1)
                & (summary_df[f"abs(Z score)"].abs() < z_threshold)
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

    # Set `Gene ID` as index
    # summary_df.set_index("Gene ID", inplace=True)

    gene_id_names_intersect = summary_df.loc[degs_intersect, "Gene ID"]
    gene_id_names_diff = summary_df.loc[degs_diff, "Gene ID"]
    gene_id_names_generic = summary_df.loc[degs_generic, "Gene ID"]

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
