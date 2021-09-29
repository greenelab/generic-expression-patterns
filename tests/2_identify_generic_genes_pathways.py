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

# # Test: Identify generic human genes on test set
#
# This notebook performs the following steps to identify generic genes
# 1. Simulates N gene expression experiments using [ponyo](https://github.com/ajlee21/ponyo)
# 2. Perform DE analysis to get association statistics for each gene
#
# In this case the DE analysis is based on the experimental design of the template experiment, described in the previous [notebook](1_process_recount2_data.ipynb). The template experiment is [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764), which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients. So the DE analysis is comparing tumor vs normal in this case.
#
# 3. For each gene, aggregate statsitics across all simulated experiments
# 4. Rank genes based on this aggregated statistic
#
# **Evaluation:**
# We want to compare our ranking using ponyo, compared to the ranking found from Crow et. al.

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

import os
import pandas as pd
import numpy as np
import pickle
from rpy2.robjects import pandas2ri
from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import process, stats, ranking

pandas2ri.activate()
np.random.seed(123)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(os.path.join(base_dir, "configs", "config_test.tsv"))

params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
NN_architecture = params["NN_architecture"]
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
count_threshold = params["count_threshold"]
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

# Load metadata file with mapping between experiments and associated samples
metadata_simulate_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "recount2_metadata.tsv"
)
metadata_delimiter = ("\t",)
experiment_id_colname = "project"

# Load pickled file
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)
# -

# ## Test: Simulation

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

# ## Test: Processing template experiment
#
# `processed_template_filename`: template data with some sample rows dropped

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

# +
# Read data
template_data = pd.read_csv(
    processed_template_filename, header=0, sep="\t", index_col=0
)

# Check samples dropped
print(template_data.shape)
assert template_data.shape[0] == 24
template_data.head()
# -

# ## Test: Processing simulation experiments

# This step modifies the following files:
# "<local_dir>/pseudo_experiments/selected_simulated_data_SRP012656_<n>.txt"
for i in range(num_runs):
    simulated_filename = os.path.join(
        local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_{i}.txt"
    )
    stats.process_samples_for_DESeq(
        simulated_filename,
        metadata_filename,
        None,
        count_threshold,
        sample_id_metadata_filename,
    )

# Check simulated files were created
sim_output1 = os.path.join(
    local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_0.txt"
)
sim_output2 = os.path.join(
    local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_1.txt"
)
assert os.path.exists(sim_output1) and os.path.exists(sim_output2)

# Check that simulated files are non-empty
assert os.path.getsize(sim_output1) > 0 and os.path.getsize(sim_output2) > 0

# +
# Check that simulated experiments are different
sim1 = pd.read_csv(sim_output1, sep="\t", index_col=0, header=0)
sim2 = pd.read_csv(sim_output2, sep="\t", index_col=0, header=0)

assert sim1.equals(sim2) == False
# -

# ## Test: Differential expression analysis

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_SRP012656_real.txt"
# get_DE_stats_DESeq(metadata_filename,
#                    project_id,
#                    processed_template_filename,
#                    "template",
#                    local_dir,
#                    "real")

# + magic_args="-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs -i base_dir" language="R"
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
#     get_DE_stats_DESeq(metadata_filename,
#                        project_id,
#                        simulated_data_filename,
#                        "simulated",
#                        local_dir,
#                        i)
# }
# -

# Check DE stats files were created
DE_output1 = os.path.join(
    local_dir, "DE_stats", "DE_stats_simulated_data_SRP012656_0.txt"
)
DE_output2 = os.path.join(
    local_dir, "DE_stats", "DE_stats_simulated_data_SRP012656_1.txt"
)
assert os.path.exists(DE_output1) and os.path.exists(DE_output2)

# Check that DE stats files are non-empty
assert os.path.getsize(DE_output1) > 0 and os.path.getsize(DE_output2) > 0

# ### Rank genes

# +
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
# -

# ### Gene summary table

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

# Some genes will have NaN's in the simulated statistics columns. These are genes that were filtered
# due to low expression and therefore the corresponding Z-score for this gene is NaN
summary_gene_ranks.isna().any()

summary_gene_ranks[summary_gene_ranks.isna().any(axis=1)]

# Create `gene_summary_filename`
gene_summary_filename = os.path.join(local_dir, "gene_summary_table.tsv")
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")

# ## Test: Compare gene ranking
# Studies have found that there are some genes that are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represents a more systematic change.
#
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# +
# Get generic genes identified by Crow et. al.
DE_prior_file = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

corr_stats, shared_ranking = ranking.compare_gene_ranking(
    summary_gene_ranks, DE_prior_file, ref_gene_col, ref_rank_col, figure_filename
)
# -

# ## Test: GSEA

# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "GSA_stats"), exist_ok=True)

# Load pathway data
hallmark_DB_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "hallmark_DB.gmt"
)

# + magic_args="-i base_dir -i template_DE_stats_filename -i hallmark_DB_filename -i statistic -o template_enriched_pathways" language="R"
# source(paste0(base_dir, '/generic_expression_patterns_modules/GSEA_analysis.R'))
#
# out_file <- paste(local_dir,
#                   "GSA_stats/GSEA_stats_template_data_",
#                   project_id,
#                   "_test.txt",
#                   sep = "")
#
# template_enriched_pathways <- find_enriched_pathways(template_DE_stats_filename, hallmark_DB_filename, statistic)
#
# template_enriched_pathways <- as.data.frame(template_enriched_pathways[1:7])
#
# write.table(template_enriched_pathways, file = out_file, row.names = F, sep = "\t")

# + magic_args="-i project_id -i local_dir -i hallmark_DB_filename -i num_runs -i statistic -i base_dir" language="R"
#
# source(paste0(base_dir,'/generic_expression_patterns_modules/GSEA_analysis.R'))
#
# # New files created: "<local_dir>/GSEA_stats/GSEA_stats_simulated_data_<project_id>_<n>.txt"
# for (i in 0:(num_runs-1)) {
#     simulated_DE_stats_file <- paste(local_dir,
#                                      "DE_stats/DE_stats_simulated_data_",
#                                      project_id,
#                                      "_",
#                                      i,
#                                      ".txt",
#                                      sep = "")
#
#     out_file <- paste(local_dir,
#                      "GSA_stats/GSEA_stats_simulated_data_",
#                      project_id,
#                      "_",
#                      i,
#                      ".txt",
#                      sep = "")
#
#     enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_filename, statistic)
#
#     # Remove column with leading edge since its causing parsing issues
#     write.table(as.data.frame(enriched_pathways[1:7]), file = out_file, row.names = F, sep = "\t")
# }
# -

# Check GSEA stats files were created
GSEA_output1 = os.path.join(
    local_dir, "GSA_stats", "GSEA_stats_simulated_data_SRP012656_0.txt"
)
GSEA_output2 = os.path.join(
    local_dir, "GSA_stats", "GSEA_stats_simulated_data_SRP012656_1.txt"
)
assert os.path.exists(DE_output1) and os.path.exists(DE_output2)

# Check that GSEA stats files are non-empty
assert os.path.getsize(GSEA_output1) > 0 and os.path.getsize(GSEA_output2) > 0

# ### Rank pathways

# +
analysis_type = "GSA"
template_GSEA_stats_filename = os.path.join(
    local_dir, "GSA_stats", f"GSEA_stats_template_data_{project_id}_test.txt"
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
# -

# ### Pathway summary table

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

# ## Test: Compare pathway ranking

# Load Powers et. al. results file
powers_rank_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "Hallmarks_qvalues_GSEAPreranked.csv"
)

# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_filename, header=0, index_col=0)
powers_rank_df.drop(["Category"], axis=1, inplace=True)

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

corr_stats = ranking.compare_pathway_ranking(
    summary_pathway_ranks, powers_rank_processed_filename, figure_filename
)
