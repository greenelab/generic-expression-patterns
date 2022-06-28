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
#     display_name: Python [conda env:generic_expression_new]
#     language: python
#     name: conda-env-generic_expression_new-py
# ---

# # Identify generic genes and pathways
#
# Studies have found that some genes are more likely to be differentially expressed even across a wide range of experimental designs. These generic genes and subsequent pathways are not necessarily specific to the biological process being studied but instead represent a more systematic change.
#
# This notebook identifies generic genes and pathways and then evaluates if those identified are consistent with published findings.
#
# **Steps to identify generic genes:**
# 1. Simulates N gene expression experiments using [ponyo](https://github.com/ajlee21/ponyo)
# 2. Perform DE analysis to get association statistics for each gene
# 3. For each gene, aggregate statsitics across all simulated experiments
# 4. Rank genes based on this aggregated statistic (i.e. log fold change, or p-value)
#
#
# **Steps to identify generic gene sets (pathways):**
# 1. Using the same simulated experiments from above, perform GSEA analysis. This analysis will determine whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list of genes, where genes are ranked by log fold change, indicating a correlation with change in expression.
# 2. For each gene set (pathway), aggregate statistics across all simulated experiments
# 3. Rank gene sets based on this aggregated statistic
#
# **Evaluation:**
# * We want to compare the ranking of genes identified using the above method with the ranking found from [Crow et. al.](https://www.pnas.org/content/pnas/116/13/6491.full.pdf), which identified a set of genes as generic based on how frequently they were found to be DE across 600 experiments
# * We want to compare the ranking of pathways identified using the above method with the ranking based on the [Powers et. al.](https://www.biorxiv.org/content/10.1101/259440v1.full.pdf) data, where ranking was determined based on the fraction of 442 experiments a pathway was found to be enriched
# * This comparison will validate our method being used as a way to automatically identify generic genes and pathways.

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2
# %matplotlib inline

import os
import sys
import pandas as pd
import numpy as np
import pickle
import glob
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
    os.path.join(base_dir, "configs", "config_human_Crow.tsv")
)

params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
NN_architecture = params["NN_architecture"]
latent_dim = params["latent_dim"]
num_runs = params["num_simulated"]
project_id = params["project_id"]
metadata_col_id = params["metadata_colname"]
mapped_template_filename = params["mapped_template_filename"]
processed_template_filename = params["processed_template_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]
scaler_filename = params["scaler_filename"]
col_to_rank_genes = params["rank_genes_by"]
col_to_rank_pathways = params["rank_pathways_by"]
statistic = params["gsea_statistic"]
logFC_name = params["DE_logFC_name"]
pvalue_name = params["DE_pvalue_name"]

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
grp_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_groups.tsv"
)

# Load metadata file with mapping between experiments and associated samples
metadata_simulate_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "experiment_sample_annotations.csv"
)
metadata_delimiter = ","
experiment_id_colname = "Experiment id"

# Load pickled file
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)

# Percentile threshold to identify generic genes
percentile_threshold = 80.0

# +
# Output files
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

pathway_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_pathway_summary_{project_id}.tsv"
)
# -

# ### Simulate experiments using selected template experiment
#
# Workflow:
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
simulate_expression_data.shift_template_experiment(
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

# ## Reverse transformation

# +
# Undo log transformation...
# CHeck that template is also in log or unlogged form
# Or log transform mapped_template instead
mapped_template = pd.read_csv(mapped_template_filename, sep="\t", index_col=0, header=0)

mapped_template_transform = (
    np.log10(mapped_template).replace(-np.inf, 0.0).replace(np.nan, 0.0)
)

mapped_template_transform.to_csv(mapped_template_filename, sep="\t")
# -

# ### Process template and simulated experiments
#
# * Remove samples not required for comparison. Since this experiment contains multiple conditions (i.e. estradiol vs EtOH at 12, 24, and 48 hrs are each considered a different comparison) being tested, we will only include those samples within the same condition.
# * Make sure ordering of samples matches metadata for proper comparison

# +
if not os.path.exists(sample_id_metadata_filename):
    sample_id_metadata_filename = None

stats.process_samples_for_limma(
    mapped_template_filename,
    grp_metadata_filename,
    processed_template_filename,
    sample_id_metadata_filename,
)

for i in range(num_runs):
    simulated_filename = os.path.join(
        local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_{i}.txt"
    )
    out_simulated_filename = os.path.join(
        local_dir,
        "pseudo_experiment",
        f"selected_simulated_data_{project_id}_{i}_processed.txt",
    )
    stats.process_samples_for_limma(
        simulated_filename,
        grp_metadata_filename,
        out_simulated_filename,
        sample_id_metadata_filename,
    )
# -

# ### Differential expression analysis
#
# The gene expression dataset is array-based so we will use Limma in this case

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i grp_metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_SRP012656_real.txt"
# get_DE_stats_limma(grp_metadata_filename,
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
    (template_DE_stats["adj.P.Val"] < 0.05) & (abs(template_DE_stats["logFC"]) > 1)
]
print(selected.shape)

# + magic_args="-i grp_metadata_filename -i project_id -i base_dir -i local_dir -i num_runs" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_SRP012656_<n>.txt"
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(local_dir,
#                                      "pseudo_experiment/selected_simulated_data_",
#                                      project_id,
#                                      "_",
#                                      i,
#                                      ".txt",
#                                      sep = "")
#
#     get_DE_stats_limma(grp_metadata_filename,
#                        project_id,
#                        simulated_data_filename,
#                        "simulated",
#                        local_dir,
#                        i)
# }
# -

# ### Rank genes

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

summary_gene_ranks.head()
# -

summary_gene_ranks.loc["CTDSP1"]

# Check if there is an NaN values, there should not be
summary_gene_ranks.isna().any()

# Create `gene_summary_fielname`
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")

# ### Compare gene ranking
# Studies have found that some genes are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represent a more systematic change.
#
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes were ranked by the number differentially expressed gene sets a gene appeared in across 600 experiments.

# +
# Get generic genes identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

corr, shared_ranking = ranking.compare_gene_ranking(
    summary_gene_ranks, DE_prior_filename, ref_gene_col, ref_rank_col, figure_filename
)
# -

# Hypergeometric test:
#
# Given $N$ number of genes with $K$ common genes in Crow et al. SOPHIE identifies $n$ genes as being common. What is the probability that $k$ of the genes identified by SOPHIE are also common in Crow et al.? What is the probability of drawing $k$ or more concordant genes?
#
# This was a way for us to quantify the correlation between SOPHIE and Crow et al common findings, since the correlation coefficient wasn't very convincing since we're considering all genes in addition to the common ones

num_Crow_genes = shared_ranking.shape[0]
num_generic_Crow_genes = shared_ranking.query(f"{ref_rank_col}>=80.0").shape[0]
num_generic_SOPHIE_genes = shared_ranking[
    shared_ranking["Percentile (simulated)"] >= percentile_threshold
].shape[0]
num_concordant_generic_genes = shared_ranking[
    (shared_ranking[ref_rank_col] >= percentile_threshold)
    & (shared_ranking["Percentile (simulated)"] >= percentile_threshold)
].shape[0]

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

# **Takeaway:**
# * Here we are comparing gene percentiles obtained from SOPHIE trained on Crow et al. vs gene percentiles obtained from manual curation using Crow et al. There is a high correlation between gene percentiles, validating the SOPHIE approach.

# ### GSEA
# **Goal:** To detect modest but coordinated changes in prespecified sets of related genes (i.e. those genes in the same pathway or share the same GO term).
#
# 1. Rank all genes using DE association statistics.
# 2. An enrichment score (ES) is defined as the maximum distance from the middle of the ranked list. Thus, the enrichment score indicates whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list (indicating a correlation with change in expression).
# 3. Estimate the statistical significance of the ES by a phenotypic-based permutation test in order to produce a null distribution for the ES (i.e. scores based on permuted phenotype)

# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "GSA_stats"), exist_ok=True)

# Load pathway data
hallmark_DB_filename = params["pathway_DB_filename"]

# + magic_args="-i base_dir -i template_DE_stats_filename -i hallmark_DB_filename -i statistic -i local_dir -o template_enriched_pathways" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/GSEA_analysis.R'))
#
# out_filename <- paste(local_dir,
#                      "GSA_stats/GSEA_stats_template_data_",
#                      project_id,
#                      "_real.txt",
#                      sep = "")
#
# template_enriched_pathways <- find_enriched_pathways(template_DE_stats_filename, hallmark_DB_filename, statistic)
# template_enriched_pathways <- as.data.frame(template_enriched_pathways[1:7])
#
# write.table(template_enriched_pathways, file = out_filename, row.names = F, sep = "\t")
# -

print(template_enriched_pathways.shape)
template_enriched_pathways[template_enriched_pathways["padj"] < 0.05].sort_values(
    by="padj"
)

# **Quick check:** Looks like enriched pathways are consistent with estradiol being estrogen hormone treatment.

# + magic_args="-i project_id -i local_dir -i hallmark_DB_filename -i num_runs -i statistic -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/GSEA_analysis.R'))
#
# # New files created: "<local_dir>/GSEA_stats/GSEA_stats_simulated_data_<project_id>_<n>.txt"
# for (i in 0:(num_runs-1)) {
#     simulated_DE_stats_filename <- paste(local_dir,
#                                      "DE_stats/DE_stats_simulated_data_",
#                                      project_id,
#                                      "_",
#                                      i,
#                                      ".txt",
#                                      sep = "")
#
#     out_filename <- paste(local_dir,
#                      "GSA_stats/GSEA_stats_simulated_data_",
#                      project_id,
#                      "_",
#                      i,
#                      ".txt",
#                      sep = "")
#
#     enriched_pathways <- find_enriched_pathways(simulated_DE_stats_filename, hallmark_DB_filename, statistic)
#
#     write.table(as.data.frame(enriched_pathways[1:7]), file = out_filename, row.names = F, sep = "\t")
# }
# -

# ### Rank pathways

analysis_type = "GSA"
template_GSEA_stats_filename = os.path.join(
    local_dir, "GSA_stats", f"GSEA_stats_template_data_{project_id}_real.txt"
)
(
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
) = ranking.process_and_rank_genes_pathways(
    template_GSEA_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_pathways,
    logFC_name,
    pvalue_name,
    "GSEA",
)

# ### Pathway summary table

# +
# Create intermediate file: "<local_dir>/gene_summary_table_<col_to_rank_pathways>.tsv"
summary_pathway_ranks = ranking.generate_summary_table(
    template_GSEA_stats_filename,
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
    col_to_rank_pathways,
    local_dir,
    "pathway",
    params,
)

summary_pathway_ranks.sort_values(by="Rank (simulated)", ascending=False).head(10)
# -

# Create `pathway_summary_filename`
summary_pathway_ranks.to_csv(pathway_summary_filename, sep="\t")

# ### Compare pathway ranking

# Studies have found that there are some pathways (gene sets) that are more likely to be significantly enriched in DEGs across a wide range of experimental designs. These generic pathways are not necessarily specific to the biological process being studied but instead represents a more systematic change.
#
# We want to compare the ability to detect these generic pathways using our method vs those found by [Powers et. al.](https://www.biorxiv.org/content/10.1101/259440v1.full.pdf) publication.  We will use the `Hallmarks_qvalues_GSEAPreranked.csv` file from https://www.synapse.org/#!Synapse:syn11806255 as a reference. The file contains the q-value (adjusted p-value) for the test: given the enrichment score (ES) of the experiment is significant compared to the null distribution of enrichment scores, where the null set is generated from permuted gene sets. For each gene set (pathway) they calculate the q-value using this test.
#
#
# To get a `reference ranking`, we calculate the fraction of experiments that a given pathway was significant (q-value <0.05) and use this rank pathways. `Our ranking` is to rank pathways based on the median q-value across the simulated experiments. We can then compare `our ranking` versus the `reference ranking.`

# Load Powers et. al. results file
powers_rank_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "Hallmarks_qvalues_GSEAPreranked.csv"
)

# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_filename, header=0, index_col=0)
powers_rank_df.drop(["Category"], axis=1, inplace=True)
print(powers_rank_df.shape)
powers_rank_df.head()

# +
# Count the number of experiments where a given pathway was found to be enriched (qvalue < 0.05)
total_num_experiments = powers_rank_df.shape[1]
frac_enriched_pathways = (powers_rank_df < 0.05).sum(axis=1) / total_num_experiments

# Rank pathways from 0-50, 50 indicating that the pathways was frequently enriched
pathway_ranks = frac_enriched_pathways.rank()

powers_rank_stats_df = pd.DataFrame(
    data={
        "Fraction enriched": frac_enriched_pathways.values,
        "Powers Rank": pathway_ranks.values,
    },
    index=powers_rank_df.index,
)
powers_rank_stats_df.sort_values(by="Powers Rank", ascending=False).head()

# +
# Save reference file for input into comparison
powers_rank_processed_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked_processed.tsv",
)

powers_rank_stats_df.to_csv(
    powers_rank_processed_filename,
    sep="\t",
)

# +
figure_filename = f"pathway_ranking_{col_to_rank_pathways}.svg"

ranking.compare_pathway_ranking(
    summary_pathway_ranks, powers_rank_processed_filename, figure_filename
)
# -

# **Takeaway:**
# * Our method ranked pathways using median adjusted p-value score across simulated experiments.
# * Powers et. al. ranked pathways based on the fraction of experiments they had adjusted p-value < 0.05.
# * Here, we compare pathway ranks obtained from SOPHIE trained on Crow et al. vs pathway ranking based on manual curation using Powers et. al datasets.
# * Despite the two datasets have a different context, we still see a good correlation.
