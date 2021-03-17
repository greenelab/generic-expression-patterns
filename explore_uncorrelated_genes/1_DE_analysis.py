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

# ## Differential expression of polyA-selected vs ribo-depleted
#
# This notebook compares 6 matched samples generated using both polyA-selected and ribo-depleted strategies. Then highlights the DEGs on the correlation plot to determine if the uncorrelated genes are those that are due to differences in protocol.

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

import os
import sys
import pandas as pd
import numpy as np

from rpy2.robjects import pandas2ri
from ponyo import utils
from generic_expression_patterns_modules import stats, ranking

pandas2ri.activate()
np.random.seed(123)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
col_to_rank_genes = params["rank_genes_by"]
project_id = params["project_id"]
# -

# ## Differential expression analysis

# **RSEM**:
# * What is RSEM doing?
# * What are the estimated count outputs?
# * gene vs transcript vs isoform level?
#
# **tximport**:
# * What is this doing?
# * Is this needed for DESeq?
#
# Note we will be rounding estimated gene counts from RSEM.
#
# Should be using tximport to get counts corrected for by gene length, but raw data isn't available
#
# https://support.bioconductor.org/p/94003/#94028
#
# **DESeq2 normalization**:
#
# Normalization is needed to perform a fair comparison -- why.
#
# The main factors considered during normalization include:
# 1. sequence depth/coverage = the number of unique reads that map to a region
# 2. gene length = longer genes tend to have more reads mapped
# 3. RNA composition = what genes are expressed and how active
#
# When performing a DE analysis, since we're comparing match gene pairs, we only need to consider sequence depth and RNA composition. DESeq2 accounts for these in its median ratio --- scale sample counts by median of ratios for sample counts vs reference pseudo count per gene
#
# How to deal with normalization and DESeq2: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#
# * Comparison
#

# +
## TO DO
## Replace integer rounding to tximport
## Update dds call to use tximport object

mapped_expression_filename = "polya_ribo_expression.tsv"
processed_expression_filename = "polya_ribo_expression_processed.tsv"
metadata_filename = "polya_ribo_sample_grouping.tsv"

stats.process_samples_for_DESeq(
    mapped_expression_filename,
    metadata_filename,
    processed_expression_filename,
)

# + magic_args="-i metadata_filename -i processed_expression_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"
# get_DE_stats_DESeq(metadata_filename,
#                    "polyA_vs_ribo",
#                    processed_expression_filename,
#                    "template",
#                    local_dir,
#                    "pbta")

# +
# Get DEGs
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", "DE_stats_template_data_polyA_vs_ribo_pbta.txt"
)

template_DE_stats = pd.read_csv(
    template_DE_stats_filename, sep="\t", header=0, index_col=0
)

selected = template_DE_stats[
    (template_DE_stats["padj"] < 0.01) & (abs(template_DE_stats["log2FoldChange"]) > 1)
]
DEGs = list(selected.index)
print(len(DEGs))
selected.head()
# -

# ## Plot

# +
# Load gene_summary_filename
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

summary_gene_ranks = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)
# -

summary_gene_ranks.loc[DEGs]

# +
# Get generic genes identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}_highlight_polyA_vs_ribo.svg"

corr, shared_ranking = ranking.compare_gene_ranking_highlight(
    summary_gene_ranks,
    DE_prior_filename,
    ref_gene_col,
    ref_rank_col,
    DEGs,
    figure_filename,
)
