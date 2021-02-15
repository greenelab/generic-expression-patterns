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

# # Coverage of eADAGE LV
#
# The goal of this notebook is to examine why genes were found to be generic. Specifically, this notebook is trying to answer the question: Are generic genes found in more multiplier latent variables compared to specific genes?
#
# The eADAGE model uses a DAE to extracts patterns of gene expression activity in the latent variables (referred to as nodes in the paper). Here we are examining the coverage of generic genes within these latent variables.
#
# **Definitions:**
# * Generic genes: Are genes that are consistently differentially expressed across multiple simulated experiments.
#
# * Other genes: These are all other non-generic genes. These genes include those that are not consistently differentially expressed across simulated experiments - i.e. the genes are specifically changed in an experiment. It could also indicate genes that are consistently unchanged (i.e. housekeeping genes)

# +
# %load_ext autoreload
# %autoreload 2

import os
import random
import textwrap
import scipy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ponyo import utils
from generic_expression_patterns_modules import lv

# +
# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "pseudomonas_analysis")

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]

project_id = params["project_id"]
# -

# Output file
nonzero_figure_filename = "nonzero_LV_coverage_eADAGE_pa.svg"
highweight_figure_filename = "highweight_LV_coverage_eADAGE_pa.svg"

# ## Load data

# Get gene summary file
summary_data_filename = os.path.join(
    data_dir, f"generic_gene_summary_{project_id}_cbrB_v_WT.tsv"
)

# +
# Load gene summary data
data = pd.read_csv(summary_data_filename, sep="\t", index_col=0, header=0)

# Check that genes are unique since we will be using them as dictionary keys below
assert data.shape[0] == len(data["Gene ID"].unique())
# -

# **Manual steps to process eADAGE data**
# 1. Data downloaded from https://zenodo.org/record/580093
# 2. Get Weight matrix (`eADAGE weight matrix.xlsx`) file
# 3. Save tab (`weight matrix`) as .csv file in `local_dir`

# Load eADAGE weight matrix
eADAGE_weight_filename = os.path.join(local_dir, "eADAGE weight matrix.csv")
eADAGE_weight = pd.read_csv(eADAGE_weight_filename, sep=",", index_col=0, header=0)

eADAGE_weight.shape

# Get a rough sense for how many genes contribute to a given LV
# (i.e. how many genes have a value > 0 for a given LV)
# Notice that eADAGE is NOT sparse
(eADAGE_weight != 0).sum().sort_values(ascending=True)

# ## Get gene data
#
# Define generic genes based on simulated gene ranking. Refer to [figure](https://github.com/greenelab/generic-expression-patterns/blob/master/pseudomonas_analysis/gene_ranking_logFC.svg) as a guide.
#
# **Definitions:**
# * Generic genes: `Percentile (simulated) >= 80`
#
# (Having a high rank indicates that these genes are consistently changed across simulated experiments.)
#
# * Other genes: `Percentile (simulated) < 80`
#
# (Having a lower rank indicates that these genes are not consistently changed across simulated experiments - i.e. the genes are specifically changed in an experiment. It could also indicate genes that are consistently unchanged.)

generic_threshold = 80
dict_genes = lv.get_generic_specific_genes(data, generic_threshold)

# +
# Check overlap between eADAGE genes and our genes
eADAGE_genes = list(eADAGE_weight.index)
our_genes = list(data.index)
shared_genes = set(our_genes).intersection(eADAGE_genes)

print(len(our_genes))
print(len(shared_genes))
# -

# Drop gene ids not used in eADAGE analysis
processed_dict_genes = lv.process_generic_specific_gene_lists(dict_genes, eADAGE_weight)

# Check numbers add up
assert len(shared_genes) == len(processed_dict_genes["generic"]) + len(
    processed_dict_genes["other"]
)

# ## Get coverage of LVs
#
# For each gene (generic or other) we want to find:
# 1. The number of LVs that gene is present
# 2. The number of LVs that the gene contributes a lot to (i.e. the gene is highly weighted within that LV)

# ### Nonzero LV coverage

dict_nonzero_coverage = lv.get_nonzero_LV_coverage(processed_dict_genes, eADAGE_weight)

# Check genes mapped correctly
assert processed_dict_genes["generic"][0] in dict_nonzero_coverage["generic"].index
assert len(dict_nonzero_coverage["generic"]) == len(processed_dict_genes["generic"])
assert len(dict_nonzero_coverage["other"]) == len(processed_dict_genes["other"])

# ### High weight LV coverage

dict_highweight_coverage = lv.get_highweight_LV_coverage_pseudomonas(
    processed_dict_genes, eADAGE_weight
)

# Check genes mapped correctly
assert processed_dict_genes["generic"][0] in dict_highweight_coverage["generic"].index
assert len(dict_highweight_coverage["generic"]) == len(processed_dict_genes["generic"])
assert len(dict_highweight_coverage["other"]) == len(processed_dict_genes["other"])

# Check high weight genes obtained are in fact at the extremes of the distribution
# Quick look at the distribution of gene weights per LV
sns.distplot(eADAGE_weight["Node2"], kde=False)
plt.yscale("log")

# ### Assemble LV coverage and plot

# +
all_coverage = []
for gene_label in dict_genes.keys():
    merged_df = pd.DataFrame(
        dict_nonzero_coverage[gene_label], columns=["nonzero LV coverage"]
    ).merge(
        pd.DataFrame(
            dict_highweight_coverage[gene_label], columns=["highweight LV coverage"]
        ),
        left_index=True,
        right_index=True,
    )
    merged_df["gene type"] = gene_label
    all_coverage.append(merged_df)

all_coverage_df = pd.concat(all_coverage)
# -

all_coverage_df = lv.assemble_coverage_df(
    processed_dict_genes, dict_nonzero_coverage, dict_highweight_coverage
)
all_coverage_df.head()

# Plot coverage distribution given list of generic coverage, specific coverage
nonzero_fig = sns.boxplot(
    data=all_coverage_df,
    x="gene type",
    y="nonzero LV coverage",
    notch=True,
    palette=["#2c7fb8", "lightgrey"],
)
nonzero_fig.set_xlabel(None)
nonzero_fig.set_xticklabels(
    ["generic genes", "other genes"], fontsize=14, fontname="Verdana"
)
nonzero_fig.set_ylabel(
    textwrap.fill("Number of LVs", width=30), fontsize=14, fontname="Verdana"
)
nonzero_fig.tick_params(labelsize=14)
nonzero_fig.set_title(
    "Number of LVs genes are present in", fontsize=16, fontname="Verdana"
)

# Notice that since our weight matrix is not sparse, all genes are present in all 300 nodes, as expected

# Plot coverage distribution given list of generic coverage, specific coverage
highweight_fig = sns.boxplot(
    data=all_coverage_df,
    x="gene type",
    y="highweight LV coverage",
    notch=True,
    palette=["#2c7fb8", "lightgrey"],
)
highweight_fig.set_xlabel(None)
highweight_fig.set_xticklabels(
    ["generic genes", "other genes"], fontsize=14, fontname="Verdana"
)
highweight_fig.set_ylabel(
    textwrap.fill("Number of LVs", width=30), fontsize=14, fontname="Verdana"
)
highweight_fig.tick_params(labelsize=14)
highweight_fig.set_title(
    "Number of LVs genes contribute highly to", fontsize=16, fontname="Verdana"
)

# ## Calculate statistics
# * Is the reduction in generic coverage significant?
# * Is the difference between generic versus other genes signficant?

# +
# Test: mean number of LVs generic genes present in vs mean number of LVs that generic gene is high weight in
# (compare two blue boxes between plots)
generic_nonzero = all_coverage_df[all_coverage_df["gene type"] == "generic"][
    "nonzero LV coverage"
].values
generic_highweight = all_coverage_df[all_coverage_df["gene type"] == "generic"][
    "highweight LV coverage"
].values

(stats, pvalue) = scipy.stats.ttest_ind(generic_nonzero, generic_highweight)
print(pvalue)

# +
# Test: mean number of LVs generic genes are high weight in vs mean number of LVs other genes high weight in
# (compare blue and grey boxes in high weight plot)
other_highweight = all_coverage_df[all_coverage_df["gene type"] == "other"][
    "highweight LV coverage"
].values
generic_highweight = all_coverage_df[all_coverage_df["gene type"] == "generic"][
    "highweight LV coverage"
].values

(stats, pvalue) = scipy.stats.ttest_ind(other_highweight, generic_highweight)
print(pvalue)

# +
# Check that coverage of other and generic genes across all LVs is NOT signficantly different
# (compare blue and grey boxes in nonzero weight plot)
other_nonzero = all_coverage_df[all_coverage_df["gene type"] == "other"][
    "nonzero LV coverage"
].values
generic_nonzero = all_coverage_df[all_coverage_df["gene type"] == "generic"][
    "nonzero LV coverage"
].values

(stats, pvalue) = scipy.stats.ttest_ind(other_nonzero, generic_nonzero)
print(pvalue)
# -

# ## Get LVs that generic genes are highly weighted in
#
# Since we are using quantiles to get high weight genes per LV, each LV has the same number of high weight genes. For each set of high weight genes, we will get the proportion of generic vs other genes. We will select the LVs that have a high proportion of generic genes to examine.

# Get proportion of generic genes per LV
prop_highweight_generic_dict = lv.get_prop_highweight_generic_genes_pseudomonas(
    processed_dict_genes, eADAGE_weight
)

proportion_generic = 0.5
generic_LV = []
for k, v in prop_highweight_generic_dict.items():
    if v > proportion_generic:
        print(k, v)
        generic_LV.append(k)

# Plot distribution of weights for these nodes
node = generic_LV[0]
lv.plot_dist_weights_pseudomonas(
    node, eADAGE_weight, shared_genes, 20, all_coverage_df, f"weight_dist_{node}.svg"
)

node = generic_LV[1]
lv.plot_dist_weights_pseudomonas(
    node, eADAGE_weight, shared_genes, 20, all_coverage_df, f"weight_dist_{node}.svg"
)

# ## Save

# +
# Save plot
nonzero_fig.figure.savefig(
    nonzero_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# Save plot
highweight_fig.figure.savefig(
    highweight_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# **Takeaway:**
# * In the first nonzero box plot, generic and other genes are present in all LVs since the eADAGE model is fully connected.
# * In the second highweight box plot, generic genes are found to contribute a lot to more LVs compared to other genes. This would indicate that generic genes are contributing alot across many LVs compared to other genes, that contribute a lot more specifically to some LVs.
