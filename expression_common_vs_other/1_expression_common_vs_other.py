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

# # Expression of common vs other genes
#
# This notebook asks if there is a technical reason for these commonly DEGs to be found as frequently differentially expressed. To answer this question, we will compare the average expression of common DEGs compared to other genes

# %load_ext autoreload
import os
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ponyo import utils
from generic_expression_patterns_modules import ranking

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
project_id = params["project_id"]
col_to_rank_genes = params["rank_genes_by"]
# -

# Read in non (0-1) normalized expression compendium
if dataset_name == "pseudomonas_analysis":
    raw_compendium_filename = params["raw_compendium_filename"]
    expression = pd.read_csv(raw_compendium_filename, sep="\t", index_col=0, header=0).T
else:
    mapped_compendium_filename = params["mapped_compendium_filename"]
    expression = pd.read_csv(
        mapped_compendium_filename, sep="\t", index_col=0, header=0
    )

expression.head()

# +
# Load gene_summary_filename
if dataset_name == "pseudomonas_analysis":
    gene_summary_filename = os.path.join(
        base_dir, dataset_name, f"generic_gene_summary_{project_id}_cbrB_v_WT.tsv"
    )
else:
    gene_summary_filename = os.path.join(
        base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
    )

summary_gene_ranks = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)
# -

summary_gene_ranks.head()

# ## Get common and other genes

# +
# Get generic genes identified by Crow et. al.
if dataset_name == "pseudomonas_analysis":
    DE_prior_filename = os.path.join(
        base_dir, dataset_name, params["reference_gene_filename"]
    )
else:
    DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}_tmp.svg"

corr, shared_ranking = ranking.compare_gene_ranking(
    summary_gene_ranks,
    DE_prior_filename,
    ref_gene_col,
    ref_rank_col,
    figure_filename,
)
# -

shared_ranking.head()

# +
# Get common gene ids
if dataset_name == "pseudomonas_analysis":
    common_ranking = shared_ranking[
        (shared_ranking["Percentile (simulated)"] > 80)
        & (shared_ranking["prop DEGs"] > 80)
    ]
    common_genes = common_ranking["gene id"]
else:
    common_ranking = shared_ranking[
        (shared_ranking["Percentile (simulated)"] > 80)
        & (shared_ranking["DE_Prior_Rank"] > 80)
    ]
    common_genes = common_ranking["Gene_Name"]


print(len(common_genes))
# -

# Get other gene ids
if dataset_name == "pseudomonas_analysis":
    other_genes = set(shared_ranking["gene id"]).difference(common_genes)
else:
    other_genes = set(shared_ranking["Gene_Name"]).difference(common_genes)
print(len(other_genes))

# ## Plot average expression

# Get average expression of genes
all_expression_mean = expression.mean()

# +
# Format df for plotting
expression_mean_toplot = pd.DataFrame(
    data={
        "all genes": np.log10(all_expression_mean),
        "common genes": np.log10(all_expression_mean[common_genes]),
        "other genes": np.log10(all_expression_mean[other_genes]),
    }
)


expression_mean_toplot.head()
# -

# Violin plot of average expression highlighing common vs other genes
f = sns.violinplot(
    data=expression_mean_toplot,
    palette=["lightgrey", "#2c7fb8", "#add8e6"],
    orient="h",
)
f.set_title(f"Average {dataset_name} expression")
f.set_xlabel("log10(average expression)")

f.get_figure().savefig(
    f"average_expression_common_vs_other_{dataset_name}.png",
    format="png",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# **Takeaway:**
#
# There doesn't appear to be a difference in expression between common genes vs other genes for this analysis. It would be interesting to continue to followup on this analysis to try to rule out other possible factors that could explain why these common DEGs are so common (for future work).
