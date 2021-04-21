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

# # Examine simulation approach
#
# **Question:** Can we separate between generic and specific genes by adding gaussian noise to simulate experiments? Does VAE approach recapitulate generic genes better than gaussian noise approach?
#
# To answer this question we will compare how well SOPHIE (VAE approach) can recapitulate manually curated generic genes (Crow et al.) compared to generic genes generated using noise approach
#
# In this notebook we will:
# 1. Generate the noise simulated experiments
# 2. Compare generic genes against Crow et al. generic genes
# 3. Compare SOPHIE vs Crow et al. results against noise vs Crow et al. results. The results for SOPHIE vs Crow et al. can be found [here](http://localhost:8888/notebooks/human_general_analysis/2_identify_generic_genes_pathways.ipynb).

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

import os
import sys
import pandas as pd
import numpy as np
import pickle
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri
from ponyo import utils
from generic_expression_patterns_modules import process, stats, ranking

pandas2ri.activate()

np.random.seed(123)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
project_id = params["project_id"]
dataset_name = params["dataset_name"]
mapped_template_filename = params["mapped_template_filename"]
processed_template_filename = params["processed_template_filename"]
num_runs = params["num_simulated"]
col_to_rank_genes = params["rank_genes_by"]
count_threshold = params["count_threshold"]
logFC_name = params["DE_logFC_name"]
pvalue_name = params["DE_pvalue_name"]

# Set mean and variance for noise distribution
mu = 0
sigma = 1000

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_groups.tsv"
)

# Percentile threshold to identify generic genes
percentile_threshold = 80.0
# -

# Output files
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_noise_model.tsv"
)

# ## Simulate data using noise approach
#
# 1. Start with template experiment
# 2. Add gaussian noise vector to each sample
# 3. Process simulated data to remove any unnecessary samples

# Create subdirectory: "<local_dir>/pseudo_experiment_noise/"
os.makedirs(os.path.join(local_dir, "pseudo_experiment_noise"), exist_ok=True)

mapped_template = pd.read_csv(mapped_template_filename, sep="\t", index_col=0, header=0)

# Simulate data by adding noise
for i in range(num_runs):
    simulated_data_filename = os.path.join(
        local_dir,
        "pseudo_experiment_noise",
        f"selected_simulated_data_{project_id}_{i}.txt",
    )

    noise = np.random.normal(mu, sigma, mapped_template.shape)

    simulated_data = mapped_template + noise

    # Set any negative counts to 0
    simulated_data[simulated_data < 0] = 0

    simulated_data.to_csv(simulated_data_filename, sep="\t")

# ### Examine distribution of template data
#
# We want to play around with the amount of noise that we add and so it would be a good idea to know what the distribution looks like for the original data

print(mapped_template.mean().mean())
sns.displot(mapped_template.mean())
plt.title("Mean gene expression for template experiment")

print(mapped_template.std().mean())
sns.displot(mapped_template.std())
plt.title("Std gene expression for template experiment")

# ## Quick check
#
# Check that we are producing distinct simulated experiments (i.e. that we are not getting the same values for each simulated experiment)

mapped_template.head()

# +
simulated_data_filename_0 = os.path.join(
    local_dir,
    "pseudo_experiment_noise",
    f"selected_simulated_data_{project_id}_0_processed.txt",
)

simulated_0 = pd.read_csv(simulated_data_filename_0, sep="\t", index_col=0, header=0)

simulated_0.head()

# +
simulated_data_filename_20 = os.path.join(
    local_dir,
    "pseudo_experiment_noise",
    f"selected_simulated_data_{project_id}_20_processed.txt",
)

simulated_20 = pd.read_csv(simulated_data_filename_20, sep="\t", index_col=0, header=0)

simulated_20.head()
# -

# ## Process template and simulated experiments
#
# * Remove samples not required for comparison
# * Make sure ordering of samples matches metadata for proper comparison
# * Make sure values are cast as integers for using DESeq
# * Filter lowly expressed genes for using DESeq

# +
if not os.path.exists(sample_id_metadata_filename):
    sample_id_metadata_filename = None

stats.process_samples_for_DESeq(
    mapped_template_filename,
    metadata_filename,
    processed_template_filename,
    count_threshold,
    sample_id_metadata_filename,
)

for i in range(num_runs):
    simulated_filename = os.path.join(
        local_dir,
        "pseudo_experiment_noise",
        f"selected_simulated_data_{project_id}_{i}.txt",
    )
    out_simulated_filename = os.path.join(
        local_dir,
        "pseudo_experiment_noise",
        f"selected_simulated_data_{project_id}_{i}_processed.txt",
    )
    stats.process_samples_for_DESeq(
        simulated_filename,
        metadata_filename,
        out_simulated_filename,
        count_threshold,
        sample_id_metadata_filename,
    )
# -

# ## Differential expression analysis
#
# The gene expression dataset is using RNA-seq so we will use DESeq2 in this case

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"
# get_DE_stats_DESeq(metadata_filename,
#                    project_id,
#                    processed_template_filename,
#                    "template",
#                    local_dir,
#                    "real")

# +
# Check number of DEGs
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)

template_DE_stats = pd.read_csv(
    template_DE_stats_filename, sep="\t", header=0, index_col=0
)

selected = template_DE_stats[
    (template_DE_stats["padj"] < 0.01) & (abs(template_DE_stats["log2FoldChange"]) > 1)
]
print(selected.shape)

# + magic_args="-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_SRP012656_<n>.txt"
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(local_dir,
#                                      "pseudo_experiment_noise/selected_simulated_data_",
#                                      project_id,
#                                      "_",
#                                      i,
#                                      "_processed.txt",
#                                      sep = "")
#
#     get_DE_stats_DESeq(metadata_filename,
#                        project_id,
#                        simulated_data_filename,
#                        "simulated",
#                        local_dir,
#                        i)
# }
# -

# ## Rank genes

analysis_type = "DE"
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

# ## Gene summary table
#
# Note: Using DESeq, genes with NaN in `Adj P-value (Real)` column are those genes flagged because of the `cooksCutoff` parameter. The cook's distance as a diagnostic to tell if a single sample has a count which has a disproportionate impact on the log fold change and p-values. These genes are flagged with an NA in the pvalue and padj columns of the result table. For more information you can read [DESeq FAQs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)

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

summary_gene_ranks.head()
# -

summary_gene_ranks.isna().any()

# Create `gene_summary_filename`
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")

# ## Compare gene ranking
# Studies have found that some genes are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represent a more systematic change.
#
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# +
# Get generic genes identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

corr, shared_ranking = ranking.compare_gene_ranking(
    summary_gene_ranks, DE_prior_filename, ref_gene_col, ref_rank_col, figure_filename
)

# +
# Hypergeometric test:
# Given N number of genes with K common genes in Crow et al.
# SOPHIE identifies n genes as being common
# What is the probability that k of the genes identified by SOPHIE
# are also common in Crow et al.? What is the probability of drawing
# k or more concordant genes?

num_Crow_genes = shared_ranking.shape[0]
num_generic_Crow_genes = shared_ranking.query(f"{ref_rank_col}>=80.0").shape[0]
num_generic_SOPHIE_genes = shared_ranking[
    shared_ranking["Percentile (simulated)"] >= percentile_threshold
].shape[0]
num_concordant_generic_genes = shared_ranking[
    (shared_ranking[ref_rank_col] >= percentile_threshold)
    & (shared_ranking["Percentile (simulated)"] >= percentile_threshold)
].shape[0]
# -

print(num_Crow_genes)
print(num_generic_Crow_genes)
print(num_generic_SOPHIE_genes)
print(num_concordant_generic_genes)

p = ss.hypergeom.sf(
    num_concordant_generic_genes,
    num_Crow_genes,
    num_generic_Crow_genes,
    num_generic_SOPHIE_genes,
)
print(p)

# **Takeaway**
# * Looks like noise and VAE can both recapitulate generic genes, which is expected.
# * Looks like template experiment already expresses generic genes (refer to other [notebook](comparisons_against_template.ipynb), so adding a small amount of noise (Normal(0,2)) will still find these generic results. This is expected, given that generic genes are "generic" because they are found across many experiments
#
# What we really want is to determine if SOPHIE can better **separate** between generic and specific genes. To do this, we would need a gold standard for what are specific genes for some experiment, which we do not have. So for now we will leave the experiment as is.
