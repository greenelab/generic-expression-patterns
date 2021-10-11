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
#     display_name: Python [conda env:generic_expression_new] *
#     language: python
#     name: conda-env-generic_expression_new-py
# ---

# # Identify generic genes and pathways
#
# This notebook allows users to find common and specific genes in their experiment of interest using their newly trained VAE model.
#
# This notebook will generate a `generic_gene_summary_<experiment id>.tsv` file that contains z-scores per gene that indicates how specific a gene is the experiment in question.

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2
# %matplotlib inline

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
# -

# ## User inputs needed
#
# User needs to define the following in the [config file](../configs/config_new_model_experiment.tsv):
#
# 1. Directory on your local machine to store intermediate and output data files generated (`local_dir`). Make sure to end with `\`.
# 2. Template experiment (`raw_template_filename`). This is the experiment you are interested in studying. This experiment is expected to be a matrix with samples as row and genes as columns (tab-delimited).
# 3. Training compendium used to train VAE, including unnormalized gene mapped version (`mapped_compendium_filename`) and normalized version (`normalized_compendium_filename`). These files should have been generated from the previous notebook.
# 4. Scaler transform (`scaler_filename`) used to normalize the training compendium. This can be found in the `data/` directory within the analysis folder.
# 5. Directory (`vae_model_dir`) containing trained VAE model (.h5 files) from the previous notebook.
# 6. Size of latent dimension (`latent_dim`).
# 7. File that maps experiment ids to the associated sample ids (`experiment_to_sample_filename`)
# 8. The delimiter used in the 'experiment_to_sample_filename' file (`metadata_delimiter`)
# 9. The column header/name that contains the experiment ids (`experiment_id_colname`)
# 10. Experiment id (`project_id`) to label newly create simulated experiments.
# 11. The column header/name in the metadatathat contains the sample ids (`sample_id_colname`)
# 12. The number of experiments to simulate (`num_simulated`)
# 13. Differential expression method to use: either 'limma' or 'desesq' (`DE_method`)
# 14. Minimum average read count to filter data by (`count_threshold`)
# 15. Name of column header (`rank_genes_by`) from DE association statistic results. This column will be use to rank genes. Select "logFC", "P.Value", "adj.P.Val", "t" if using Limma. Select "log2FoldChange", "pvalue", "padj" if using DESeq.
# 16. `DE_logFC_name` is either "logFC" (Limma) or "log2FoldChange" (DESeq). This is used for plotting volcano plots.
# 17. `DE_pvalue_name` is either "adj.P.Val" (Limma) or "padj" (DESeq). This is used for plotting volcano plots.
#
# The remaining parameters within the `config` file specify filenames that are intermediate data files that will be generated when SOPHIE runs.
#
# The user also needs to provide metadata files:
#
# These files should be located in `data/metadata/` directory.
# 1. `<experiment id>_process_samples.tsv` contains 2 columns tab-delimited (sample ids, label that indicates if the sample is kept or removed). See [example](data/metadata/cis-gem-par-KU1919_process_samples.tsv). **Note: This file is not required if the user wishes to use all the samples in the template experiment file.**
# 2. `<experiment id>_groups.tsv` contains 2 columns tab-delimited: sample ids, group label to perform DE analysis. See [example](data/metadata/cis-gem-par-KU1919_groups.tsv)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_new_model_experiment.tsv")
)
params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
NN_architecture = params["NN_architecture"]
latent_dim = params["latent_dim"]
num_runs = params["num_simulated"]
metadata_simulate_filename = params["experiment_to_sample_filename"]
metadata_delimiter = params["metadata_delimiter"]
experiment_id_colname = params["experiment_id_colname"]
project_id = params["project_id"]
metadata_col_id = params["sample_id_colname"]
raw_template_filename = params["raw_template_filename"]
processed_template_filename = params["processed_template_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]
scaler_filename = params["scaler_filename"]
method = params["DE_method"]
col_to_rank_genes = params["rank_genes_by"]
logFC_name = params["DE_logFC_name"]
pvalue_name = params["DE_pvalue_name"]

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
# -

# Output files
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

# ## Simulate experiments using selected template experiment
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

# +
simulated_filename = os.path.join(
    local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_1.txt"
)

test = pd.read_csv(simulated_filename, sep="\t", index_col=0, header=0)
# -

test.head()

# ## Process template and simulated data
#
# * Remove samples not required for comparison.
# * Make sure ordering of samples matches metadata for proper comparison

# +
if not os.path.exists(sample_id_metadata_filename):
    sample_id_metadata_filename = None

if method == "deseq":
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
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt",
        )
        out_simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}_processed.txt",
        )
        stats.process_samples_for_DESeq(
            simulated_filename,
            metadata_filename,
            out_simulated_filename,
            count_threshold,
            sample_id_metadata_filename,
        )
else:
    stats.process_samples_for_limma(
        mapped_template_filename,
        metadata_filename,
        processed_template_filename,
        sample_id_metadata_filename,
    )

    for i in range(num_runs):
        simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt",
        )
        stats.process_samples_for_limma(
            simulated_filename,
            metadata_filename,
            None,
            sample_id_metadata_filename,
        )
# -

# ## Differential expression analysis

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir -i method" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"
# if (method == "deseq"){
#     get_DE_stats_DESeq(
#         metadata_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }
# else{
#     get_DE_stats_limma(
#         metadata_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }

# + magic_args="-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs -i method" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_<project_id>_<n>.txt"
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(
#         local_dir,
#         "pseudo_experiment/selected_simulated_data_",
#         project_id,
#         "_",
#         i,
#         "_processed.txt",
#         sep = ""
#     )
#     if (method == "deseq"){
#         get_DE_stats_DESeq(
#             metadata_filename,
#             project_id,
#             simulated_data_filename,
#             "simulated",
#             local_dir,
#             i
#             )
#     }
#     else {
#         get_DE_stats_limma(
#             metadata_filename,
#             project_id,
#             simulated_data_filename,
#             "simulated",
#             local_dir,
#             i
#             )
#         }
#     }
# -

# ## Rank genes
# Genes are ranked by their "generic-ness" - how frequently these genes are changed across the simulated experiments using user-specific test statistic (i.e. log2 fold change).

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

# ## Summary table
#
# * Gene ID: Gene identifier (hgnc symbols for human data or PA number for *P. aeruginosa* data)
# * (Real): Statistics for template experiment
# * (Simulated): Statistics across simulated experiments
# * Number of experiments: Number of simulated experiments
# * Z-score: High z-score indicates that gene is more changed in template compared to the null set of simulated experiments (high z-score = highly specific to template experiment)
# * Percentile (simulated): percentile rank of the median(abs(log fold change)). So its the median absolute change for that gene across the 25 simulated experiments that is then converted to a percentile rank from 0 - 100. Where a higher percentile indicates that the gene was highly changed frequently and would suggest that the gene is more commonly DE.
# * Percent DE (simulated): the fraction of the simulated experiments in which that gene was found to be DE using (log fold change > 1 and adjusted p-value < 0.05). _Note:_ you may find that many genes have a 0 fraction. This is because there is some compression that happens when pushing data through the VAE so the variance of the simulated experiments is lower compared to the real experiment. We are aware of this limitation in the VAE and are looking at how to improve the variance and biological signal captured by the VAE, however we were still able to demonstrate that for now the VAE is able to simulate realistic looking biological experiments in our previous [paper](https://academic.oup.com/gigascience/article/9/11/giaa117/5952607).
#
#
# **Note:**
# * If using DESeq, genes with NaN in only the `Adj P-value (Real)` column are those genes flagged because of the `cooksCutoff` parameter. The cook's distance as a diagnostic to tell if a single sample has a count which has a disproportionate impact on the log fold change and p-values. These genes are flagged with an NA in the pvalue and padj columns of the result table.
#
# * If using DESeq with count threshold, some genes may not be present in all simulated experiments (i.e. the `Number of experiments (simulated)` will not equal the number of simulated experiments you specified in the beginning. This pre-filtering will lead to some genes found in few simulated experiments and so the background/null set for that gene is not robust. Thus, the user should sort by both z-score and number of experiments to identify specific expressed genes.
#
# * If using DESeq without count threshold, some genes may still not be present in all simulated experiments (i.e. the `Number of experiments (simulated)`  will not equal the number of simulated experiments you specified in the beginning. If the gene is 0 expressed across all samples and thus automatically given an NA in `log fold change, adjusted p-value` columns. Thus, the user should sort by both z-score and number of experiments to identify specific expressed genes.
#
# For more information you can read [DESeq FAQs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)

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

# Check if there is an NaN values, there should not be
summary_gene_ranks.isna().any()

# Create `gene_summary_filename`
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")
