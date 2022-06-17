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

# # SOPHIE vs traditional DE
#
# The goal of this experiment is to determine how often SOPHIE ranks specific genes over generic genes compared to using traditional DE analysis

# +
# %load_ext autoreload
# %autoreload 2
import os
import pickle
import numpy as np
import pandas as pd
import random
from ponyo import utils

np.random.seed(1)
# -

# ## Create simulated data
#
# 1. Each perturbation experiment has 4 samples (2 perturbed vs 2 control) with 1000 genes. Create initial expression profiles for 4 samples by drawing from a gaussian with a mean/sd for each gene.
# 2. Say that there are 100 generic genes. The generic genes will have the same scalar value added for the “perturbed” samples
# 3. Each perturbation experiment will have 10 “specific” genes perturbed in addition to the “generic” ones. Select specific genes randomly . Then the specific genes will have the same scalar value added for the “perturbed” samples
# 4. Repeat this process 90 times to get 90 experiments

# +
# User params
num_genes = 1000
num_samples = 4
num_generic_genes = 100
num_specific_genes = 10
num_experiments = 90

# Range of counts
min_expression = 0
max_expression = 1000

# params for distribution for scaler used to shift generic and specific genes
mean_scaler = 0
var_scaler = 10

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = "config_sophie_vs_trad.tsv"

params = utils.read_config(config_filename)

# +
# Load config params

# Local directory to store intermediate files
local_dir = params["local_dir"]

# Un-normalized compendium filename
raw_compendium_filename = params["raw_compendium_filename"]
# -

# Pickle files to save
generic_gene_ids_filename = "generic_gene_ids.pickle"

# +
# Setup variables

# Sample ids (2 controls, 2 perturbed)
sample_ids = ["0_control", "1_control", "2_perturb", "3_perturb"]

# Make range of numbers into string to use as gene ids
gene_ids = [f"G_{i}" for i in range(num_genes)]

# Randomly select generic gene ids
generic_gene_ids = random.sample(gene_ids, num_generic_genes)
print(generic_gene_ids)


# -

# ## Supporting functions

# Function to make an individual experiment
def run_make_experiment(
    num_samples,
    num_genes,
    min_expression,
    max_expression,
    sample_ids,
    all_gene_ids,
    generic_gene_ids,
    specific_gene_ids,
    generic_scaler,
    specific_scaler,
):
    experiment_data = {}
    for sample in range(num_samples):

        # sample_profile = np.zeros(num_genes)
        # for j in range(num_genes):
        #    sample_profile[j] = random.randint(min_expression, max_expression)

        # Randomly select a different probability for each sample
        p = random.uniform(0.0, 1.0)
        sample_profile = np.random.negative_binomial(100, p, num_genes)

        # Create dictionary to define dataframe
        experiment_data[sample_ids[sample]] = sample_profile

    # Create experiment dataframe
    experiment_df = pd.DataFrame(experiment_data).T

    # Set column header
    experiment_df.columns = all_gene_ids

    # Perturb generic genes by scaler
    # Randomly select a scaler
    # Only add scaler to perturbed samples
    experiment_df.loc[
        experiment_df.index.str.contains("perturb"), generic_gene_ids
    ] += generic_scaler

    # Perturb specific genes by a different scaler
    experiment_df.loc[
        experiment_df.index.str.contains("perturb"), specific_gene_ids
    ] += specific_scaler

    return experiment_df


# Make multiple experiments
def make_experiments(
    num_experiments,
    num_samples,
    num_genes,
    mean_expression,
    var_expression,
    sample_ids,
    all_gene_ids,
    generic_gene_ids,
    num_specific_genes,
):

    expression_df = pd.DataFrame()
    specific_gene_id_lst = []

    for i in range(num_experiments):

        # Randomly select specific genes from the pool of (900)remaining genes
        # Select without replacement
        remaining_gene_ids = list(set(all_gene_ids).difference(generic_gene_ids))
        specific_gene_ids = random.sample(remaining_gene_ids, num_specific_genes)

        # Save specific gene ids for reference later
        specific_gene_id_lst.append(specific_gene_ids)

        # Randomly select scaler for generic and specific genes from the same distribution
        scale_factors = np.random.randint(mean_scaler, var_scaler, size=2)
        generic_scaler = scale_factors[0]
        specific_scaler = scale_factors[1]

        experiment_df = run_make_experiment(
            num_samples,
            num_genes,
            mean_expression,
            var_expression,
            sample_ids,
            all_gene_ids,
            generic_gene_ids,
            specific_gene_ids,
            generic_scaler,
            specific_scaler,
        )

        # Concatenate experiments
        expression_df = pd.concat([expression_df, experiment_df])

    # Try to reset index to see if this makes a difference
    # NOTE: VAE don't train when sample indices are identical, not sure why
    if num_experiments > 1:
        expression_df = expression_df.reset_index(drop=True)

    return expression_df, specific_gene_id_lst


# ## Make a template experiment

for i in range(10):
    template_experiment, template_specific_gene_ids = make_experiments(
        1,
        num_samples,
        num_genes,
        min_expression,
        max_expression,
        sample_ids,
        gene_ids,
        generic_gene_ids,
        num_specific_genes,
    )

    # Save template experiment
    raw_template_filename = f"/home/alexandra/Documents/Data/Generic_expression_patterns/reviewer_experiment/raw_template_{i}.tsv"
    template_experiment.to_csv(raw_template_filename, sep="\t")

    # Pickle specific gene ids
    template_specific_gene_ids_filename = f"/home/alexandra/Documents/Data/Generic_expression_patterns/reviewer_experiment/template_specific_gene_ids_{i}.pickle"
    with open(template_specific_gene_ids_filename, "wb") as pkl_fh:
        pickle.dump(template_specific_gene_ids[0], pkl_fh, protocol=3)

print(template_experiment.shape)
template_experiment.head()

template_experiment[generic_gene_ids]

template_specific_gene_ids[0]

# ## Make compendium

compendium, compendium_specific_ids = make_experiments(
    num_experiments,
    num_samples,
    num_genes,
    min_expression,
    max_expression,
    sample_ids,
    gene_ids,
    generic_gene_ids,
    num_specific_genes,
)

print(compendium.shape)
compendium.head()

compendium_specific_ids

# +
# Save
compendium.to_csv(raw_compendium_filename, sep="\t")

# Save generic genes
# Pickle `scaler` as `scaler_filename` on disk
with open(generic_gene_ids_filename, "wb") as pkl_fh:
    pickle.dump(generic_gene_ids, pkl_fh, protocol=3)
