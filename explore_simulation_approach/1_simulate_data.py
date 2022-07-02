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
import seaborn as sns
from ponyo import utils

np.random.seed(2)

# +
## Need to make sure there is a consistent difference in all experiments in the compendium
# -

# ## Create simulated data
#
# 1. Each perturbation experiment has 8 samples (4 perturbed vs 4 control) with 1000 genes. Create initial expression profiles for 4 samples by drawing from a gaussian with a mean/sd for each gene.
# 2. Say that there are 100 generic genes selected from the 1000 genes. The generic genes will have the some scalar value added for the “perturbed” samples
# 3. Each perturbation experiment will have 10 “specific” genes are randomly selected from the remaining 900 non-generic genes. Then the specific genes will have the some scalar value added for the “perturbed” samples
# 4. Repeat this process 90 times to get 90 experiments

# User params
num_genes = 1000
num_samples = 8
num_generic_genes = 100
num_specific_genes = 10
num_experiments = 90

"""r = random.randint(5,50)
print(r)
p = random.uniform(0.0, 1.0)
print(p)
backgrd_dist = np.random.negative_binomial(27, p, 1000)
#pp = random.uniform(0.0, 0.2)
#print(pp)
perturb_dist = np.random.negative_binomial(20, p, 1000)"""

"""sns.distplot(backgrd_dist, kde=False)
sns.distplot(perturb_dist, kde=False)"""

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

# Sample ids (4 controls, 4 perturbed)
sample_ids = [
    "0_control",
    "1_control",
    "2_control",
    "3_control",
    "4_perturb",
    "5_perturb",
    "6_perturb",
    "7_perturb",
]

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
    sample_ids,
    all_gene_ids,
    generic_gene_ids,
    specific_gene_ids,
    r,
):
    experiment_data = {}
    for gene in range(num_genes):

        # Simulate expression by sampling from a Negative Binomial distribution
        # NB(r,p)
        # r: number of successes until we stop trials
        # p: success rate, the probability that each trial is a success
        # Within an experiment, each trial can have 2 outcomes = success or failure
        # Within a trial, the probability of success is p
        # NB()
        # Randomly select a different success probability for each gene so that
        # each gene has a different rate at which its expressed
        p = random.uniform(0.0, 1.0)
        gene_profile = np.random.negative_binomial(r, p, num_samples)

        # Create dictionary to define dataframe
        experiment_data[f"G_{gene}"] = gene_profile

    # Create experiment dataframe
    experiment_df = pd.DataFrame(experiment_data)

    # Set index
    experiment_df.index = sample_ids

    return experiment_df


# Make multiple experiments
def make_experiments(
    num_experiments,
    num_samples,
    num_genes,
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

        # Set a different r per experiment to try to add more
        # variance to the compendium
        # Issue: the VAE is compressing the DE signal in the template
        # experiment too much so the simulated experiments don't retain this
        # DE signal at all.
        # Hypothesis: Generating more variance in the compendium experiments
        # may reduce the compression
        r = random.randint(5, 50)

        experiment_df = run_make_experiment(
            num_samples,
            num_genes,
            sample_ids,
            all_gene_ids,
            generic_gene_ids,
            specific_gene_ids,
            r,
        )

        # Randomly select scaler for generic and specific genes to be the median
        # of the expression values so that the signal is strong

        # p = random.uniform(0.0, 0.2)
        # generic_scaler = np.random.negative_binomial(100, p)
        # specific_scaler = np.random.negative_binomial(50, p)
        generic_scaler = experiment_df.median().quantile(0.98)
        specific_scaler = generic_scaler
        print("generic scaler", generic_scaler)
        print("specific scaler", specific_scaler)

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

# ## Make compendium

compendium, compendium_specific_ids = make_experiments(
    num_experiments,
    num_samples,
    num_genes,
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
with open(generic_gene_ids_filename, "wb") as pkl_fh:
    pickle.dump(generic_gene_ids, pkl_fh, protocol=3)
