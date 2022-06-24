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

# # Compare generic genes
#
# The goal of this notebook is to compare the generic genes found using the same template experiment run two times and 2 different recount2 template experiments.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline

import os
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from ponyo import utils

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]

project_id1 = "SRP012656"
project_id2 = "SRP061689"

# +
# Get data directory containing gene summary data
data_dir = os.path.join(base_dir, "human_general_analysis")

# Get gene ranking files
gene_ranking_filename1 = os.path.join(
    data_dir, f"generic_gene_summary_{project_id1}.tsv"
)
gene_ranking_filename1_run2 = os.path.join(
    data_dir, f"generic_gene_summary_{project_id1}_run2.tsv"
)
gene_ranking_filename2 = os.path.join(
    data_dir, f"generic_gene_summary_{project_id2}.tsv"
)

# Get template data
template_filename1 = os.path.join(
    data_dir, "data", f"processed_recount2_template_{project_id1}.tsv"
)
template_filename2 = os.path.join(
    data_dir, "data", f"processed_recount2_template_{project_id2}.tsv"
)
# -

# ## Correlation between rankings between same experiment
#
# Here we compare gene ranking after running SOPHIE 2 times using the same template experiment but different seeds.

# Load gene ranking
gene_ranking_summary1 = pd.read_csv(
    gene_ranking_filename1, sep="\t", index_col=0, header=0
)
gene_ranking_summary1_run2 = pd.read_csv(
    gene_ranking_filename1_run2, sep="\t", index_col=0, header=0
)

# Get simulated ranking
gene_ranking1 = (
    gene_ranking_summary1["Rank (simulated)"].rename("Rank 1").to_frame("Rank 1")
)
gene_ranking1_run2 = (
    gene_ranking_summary1_run2["Rank (simulated)"].rename("Rank 2").to_frame("Rank 2")
)

# +
# Scale ranking to percentile (0,100)
scaler = MinMaxScaler(feature_range=(0, 100))

gene_ranking1["Percentile 1"] = scaler.fit_transform(
    np.array(gene_ranking1["Rank 1"]).reshape(-1, 1)
)

gene_ranking1_run2["Percentile 2"] = scaler.fit_transform(
    np.array(gene_ranking1_run2["Rank 2"]).reshape(-1, 1)
)

gene_ranking1_run2.head()
# -

# Combine ranking
gene_ranking_same_combined = pd.concat(
    [gene_ranking1["Percentile 1"], gene_ranking1_run2["Percentile 2"]], axis=1
)

print(gene_ranking_same_combined.shape)
gene_ranking_same_combined.head()

# Check for NAs
gene_ranking_same_combined[pd.isnull(gene_ranking_same_combined).any(axis=1)]

# +
# Plot correlation between ranking
r, p = stats.spearmanr(
    gene_ranking_same_combined["Percentile 1"],
    gene_ranking_same_combined["Percentile 2"],
)
print(r, p)

fig = sns.jointplot(
    data=gene_ranking_same_combined,
    x="Percentile 1",
    y="Percentile 2",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

# Make axis thicker
for _, s in fig.fig.axes[0].spines.items():
    s.set_linewidth(1.5)

# Remove marginal axes
fig.fig.axes[1].set_visible(False)
fig.fig.axes[2].set_visible(False)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

fig.set_axis_labels(
    f"Percentile in {project_id1} seed 1",
    f"Percentile in {project_id1} seed 2",
    fontsize=20,
    fontname="Verdana",
)
cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes", fontsize=12)
cb.ax.tick_params(labelsize=12)

output_figure_filename = "concordance_between_same_recount2_templates.svg"
fig.savefig(
    output_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# **Takeaway:**
# * Running SOPHIE twice using the same template experiment will generate 2 different sets of simulated experiments.
# * Since the template experiment is the same, these 2 sets of simulated experiments will have the same experimental design structure/biological context
# * As expected, the concordance is very high especially for high ranked and low ranked genes. The genes in the middle rank are more sensitive to changes so you don't get as clear of a signal compared to the extreme ranked genes.

# ## Correlation between rankings between 2 different experiments
#
# Here we compare gene ranking generated by SOPHIE using 2 different template experiments.

# Load gene ranking
gene_ranking_summary2 = pd.read_csv(
    gene_ranking_filename2, sep="\t", index_col=0, header=0
)

# Get simulated ranking
gene_ranking1 = (
    gene_ranking_summary1["Rank (simulated)"].rename("Rank 1").to_frame("Rank 1")
)
gene_ranking2 = (
    gene_ranking_summary2["Rank (simulated)"].rename("Rank 2").to_frame("Rank 2")
)

# +
# Scale ranking to percentile (0,100)
scaler = MinMaxScaler(feature_range=(0, 100))

gene_ranking1["Percentile 1"] = scaler.fit_transform(
    np.array(gene_ranking1["Rank 1"]).reshape(-1, 1)
)

gene_ranking2["Percentile 2"] = scaler.fit_transform(
    np.array(gene_ranking2["Rank 2"]).reshape(-1, 1)
)

gene_ranking2.head()
# -

# Combine ranking
gene_ranking_diff_combined = pd.concat(
    [gene_ranking1["Percentile 1"], gene_ranking2["Percentile 2"]], axis=1
)

print(gene_ranking_diff_combined.shape)
gene_ranking_diff_combined.head()

# Check for NAs
gene_ranking_diff_combined[pd.isnull(gene_ranking_diff_combined).any(axis=1)]

# +
# Plot correlation between ranking
r, p = stats.spearmanr(
    gene_ranking_diff_combined["Percentile 1"],
    gene_ranking_diff_combined["Percentile 2"],
)
print(r, p)

fig = sns.jointplot(
    data=gene_ranking_diff_combined,
    x="Percentile 1",
    y="Percentile 2",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

# Make axis thicker
for _, s in fig.fig.axes[0].spines.items():
    s.set_linewidth(1.5)

# Remove marginal axes
fig.fig.axes[1].set_visible(False)
fig.fig.axes[2].set_visible(False)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

fig.set_axis_labels(
    f"Percentile in {project_id1}",
    f"Percentile in {project_id2}",
    fontsize=20,
    fontname="Verdana",
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes", fontsize=12)
cb.ax.tick_params(labelsize=12)

output_figure_filename = "concordance_between_diff_recount2_templates.svg"
fig.savefig(
    output_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# **Takeaway:**
#
# * Looks like there is good concordance between highly ranked genes (i.e. generic genes)
# * By comparison if we run SOPHIE using two different template experiments, there are genes in the off-diagonal regions that might indicate that there are generic within the given context of the specific experiment.
# * In general, the genes in the middle rank are more sensitive to changes so you don't get as clear of a signal compared to the highest rank genes.

# ## Examine gene expression data

# Read expression data
template_1 = pd.read_csv(template_filename1, sep="\t", index_col=0, header=0)
template_2 = pd.read_csv(template_filename2, sep="\t", index_col=0, header=0)

# +
# Get concordance genes
concordant_genes = list(
    gene_ranking_diff_combined[
        (gene_ranking_diff_combined["Percentile 1"] > 80)
        & (gene_ranking_diff_combined["Percentile 2"] > 80)
    ].index
)

# Get disconcordant genes
discordant_genes = set(gene_ranking_diff_combined.index).difference(concordant_genes)

# +
# Distribution of concordant genes in template experiment 1
template1_mean = template_1.mean()

print(
    "Percent concordant genes with 0 expression in template 1:",
    len(template1_mean[concordant_genes].loc[template1_mean[concordant_genes] == 0])
    / len(template1_mean[concordant_genes]),
)

print(
    "Percent nonzero concordant genes in template 1:",
    len(
        template1_mean[concordant_genes].loc[
            (template1_mean[concordant_genes] > 0)
            & (template1_mean[concordant_genes] < 1000)
        ]
    )
    / len(template1_mean[concordant_genes]),
)

f1 = sns.distplot(template_1.mean()[concordant_genes], kde=False)
f1.set_title(f"Expression of concordant genes in {project_id1}")
f1.set_xlabel("log(gene expression)")
f1.set_ylabel("log(count)")
f1.set(xscale="log", yscale="log")

# +
# Distribution of concordant genes in template experiment 2
template2_mean = template_2.mean()
print(
    "Percent concordant genes with 0 expression in template 2:",
    len(template2_mean[concordant_genes].loc[template2_mean[concordant_genes] == 0])
    / len(template2_mean[concordant_genes]),
)

print(
    "Percent nonzero concordant genes in template 2:",
    len(
        template2_mean[concordant_genes].loc[
            (template2_mean[concordant_genes] > 0)
            & (template2_mean[concordant_genes] < 1000)
        ]
    )
    / len(template2_mean[concordant_genes]),
)

# There are more 0 expressed genes in this template experiment
f2 = sns.distplot(template_2.mean()[concordant_genes], kde=False)
f2.set_title(f"Expression of concordant genes in {project_id2}")
f2.set_xlabel("log(gene expression)")
f2.set_ylabel("log(count)")
f2.set(xscale="log", yscale="log")

# +
# Distribution of discordant gense in template experiment 1
template1_mean = template_1.mean()

print(
    "Percent discordant genes with 0 expression in template 1:",
    len(template1_mean[discordant_genes].loc[template1_mean[discordant_genes] == 0])
    / len(template1_mean[discordant_genes]),
)

print(
    "Percent nonzero discordant genes in template 1:",
    len(
        template1_mean[discordant_genes].loc[
            (template1_mean[discordant_genes] > 0)
            & (template1_mean[discordant_genes] < 1000)
        ]
    )
    / len(template1_mean[discordant_genes]),
)

print(
    len(template1_mean[discordant_genes].loc[template1_mean[discordant_genes] > 0])
    / len(template1_mean[discordant_genes])
)
f3 = sns.distplot(template_1.mean()[discordant_genes], kde=False)
f3.set_title(f"Expression of discordant genes in {project_id1}")
f3.set_xlabel("log(gene expression)")
f3.set_ylabel("log(count)")
f3.set(xscale="log", yscale="log")

# +
# Distribution of discordant genes in template experiment 2
template2_mean = template_2.mean()

print(
    "Percent discordant genes with 0 expression in template 2:",
    len(template2_mean[discordant_genes].loc[template2_mean[discordant_genes] == 0])
    / len(template2_mean[discordant_genes]),
)

print(
    "Percent nonzero discordant genes in template 2:",
    len(
        template2_mean[discordant_genes].loc[
            (template2_mean[discordant_genes] > 0)
            & (template2_mean[discordant_genes] < 1000)
        ]
    )
    / len(template2_mean[discordant_genes]),
)

f4 = sns.distplot(template_2.mean()[discordant_genes], kde=False)
f4.set_title(f"Expression of discordant genes in {project_id2}")
f4.set_xlabel("log(gene expression)")
f4.set_ylabel("log(count)")
f4.set(xscale="log", yscale="log")
# -

# **Takeaway:**
#
# Doesn't appear to be much of a difference between the distribution of average gene expression values for these two experiments.
#
# Theoretically, I would expect the scenario where a gene is lowly expressed in the context of template experiment 1 and therefore not found to be generic. But this same gene could be found to be generic in the context of template experiment 2 if it is more expressed. Its possible that differences in gene expression distribution can change which genes are found to be generic given that the simulation is producing experiments with a similar context.
#
# In this case, despite having similar gene expression distributions there are still many differences in gene ranking. This suggests to me that level of gene expression activity doesn't matter as much as the overall patterns perhaps.
#
# Overall we observe a slight shift showing that concordant genes are more lowly expressed compared to discordant genes, but most genes are still predominantly lowly gene expression. If most genes have expression levels very close to 0, then small fluctuations in the expression of some genes could lead to large changes in rank without changing the overall expression distribution.
