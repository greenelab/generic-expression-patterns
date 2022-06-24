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

# # Plot for E-GEOD-33245 (cbrB vs WT)
#
# This notebook is creating volcano plots that highlight the ArgR regulon, which are the genes selected based on the SOPHIE analysis to be specific. These plots will correspond to the wet lab experiments conducted by the Hogan Lab.

# +
# %load_ext autoreload
# %matplotlib inline

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ponyo import utils

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)

params = utils.read_config(config_filename)
# -

# Load config params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
project_id = params["project_id"]
num_simulated = params["num_simulated"]

# Load DE stats
template_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_cbrB_v_WT.tsv"
)

# Load argR genes
# These were the selected specific genes used in the wet lab experiments
argR_genes_filename = os.path.join("data", "ArgR.csv")

argR_genes = list(pd.read_csv(argR_genes_filename, index_col=0, header=0).index)


# ## Plotting functions

def make_volcano_template_highlight_genelist(
    template_summary_filename,
    project_id,
    y_stat,
    genes_to_highlight,
    output_figure_filename,
):
    """
    This function creates volcano plot of template experiment
    highlighting traditional DEGs

    Arguments
    ----------
    template_DE_stats_filename: str
        File containing summary statistics for template experiment
    project_id: str
        Experiment identifier
    y_stat: 'adj.P.Val' or 'Z score'
    genes_to_highlight: list
        List of genes ids to highlight
    output_figure_filename: str
        File to save figure to
    """

    # Read template DE stats
    template_summary_df = pd.read_csv(
        template_summary_filename, sep="\t", index_col=0, header=0
    )

    if y_stat == "adj.P.Val":
        # Take -log10 of adjusted p-value
        template_summary_df["padj_log10"] = -np.log10(
            template_summary_df["Adj P-value (Real)"]
        )

    # Label genes in input list
    template_summary_df["gene group"] = "none"
    template_summary_df.loc[genes_to_highlight, "gene group"] = "ArgR genes"

    # Plot
    # Note: Tried to use
    # hue="gene group",
    # hue_order=["none", "ArgR genes"]
    # parameters but
    # the genes I want to highlight still seem to be covered
    # for some reason so I have switched to plot them separately instead
    if y_stat == "adj.P.Val":
        f = sns.scatterplot(
            data=template_summary_df[template_summary_df["gene group"] == "none"],
            x="logFC (Real)",
            y="padj_log10",
            hue="gene group",
            alpha=0.5,
            palette=["lightgrey"],
            linewidth=0,
            legend=False,
        )
        f = sns.scatterplot(
            data=template_summary_df[template_summary_df["gene group"] == "ArgR genes"],
            x="logFC (Real)",
            y="padj_log10",
            hue="gene group",
            alpha=0.5,
            palette=["red"],
            linewidth=0,
        )

        # Add traditional thresholds
        f.axhline(-np.log10(0.05), c="black", lw=0.7, ls="--")
        f.axvline(1, c="black", lw=0.7, ls="--")
        f.axvline(-1, c="black", lw=0.7, ls="--")

        # Move location of legend
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

        f.set_ylabel(
            r"-log$_{10}$ (FDR adjusted p-value)", fontsize=14, fontname="Verdana"
        )
    else:
        f = sns.scatterplot(
            data=template_summary_df[template_summary_df["gene group"] == "none"],
            x="logFC (Real)",
            y="Z score",
            hue="gene group",
            alpha=0.5,
            palette=["lightgrey"],
            linewidth=0,
            legend=False,
        )
        f = sns.scatterplot(
            data=template_summary_df[template_summary_df["gene group"] == "ArgR genes"],
            x="logFC (Real)",
            y="Z score",
            hue="gene group",
            alpha=0.5,
            palette=["red"],
            linewidth=0,
        )
        f.set_ylabel("Specificity score (z-score)", fontsize=14, fontname="Verdana")

        # Move location of legend
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

    f.set_xlabel(r"log$_2$ Fold Change", fontsize=14, fontname="Verdana")
    f.set_title(f"Template experiment ({project_id})", fontsize=16, fontname="Verdana")

    f.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


len(argR_genes)

# ## logFC vs p-value

make_volcano_template_highlight_genelist(
    template_summary_filename,
    project_id,
    "adj.P.Val",
    argR_genes,
    os.path.join(local_dir, f"template_traditional_volcano_ArgR_{project_id}.svg"),
)

# ## logFC vs z-score

make_volcano_template_highlight_genelist(
    template_summary_filename,
    project_id,
    "z-score",
    argR_genes,
    os.path.join(local_dir, f"template_zscore_volcano_ArgR_{project_id}.svg"),
)

# ## logFC vs z-score difference

# Load summary tables for cbrB vs WT and crc vs WT
cbrB_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_cbrB_v_WT.tsv"
)
crc_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_crc_v_WT.tsv"
)

# Read summary tables
cbrB_summary_df = pd.read_csv(cbrB_summary_filename, sep="\t", index_col=0, header=0)
crc_summary_df = pd.read_csv(crc_summary_filename, sep="\t", index_col=0, header=0)

cbrB_summary_df.head()

crc_summary_df.head()

# Select logFC(Real) and z score columns
cbrB_select = cbrB_summary_df[["Gene Name", "logFC (Real)", "Z score"]]
crc_select = crc_summary_df[["logFC (Real)", "Z score"]]

print(cbrB_select.shape)
print(crc_select.shape)

# Merge on gene id
cbrB_crc_df = cbrB_select.merge(
    crc_select, left_index=True, right_index=True, suffixes=["_cbrB", "_crc"]
)
print(cbrB_crc_df.shape)
cbrB_crc_df.head()

# Calculate the difference in z score
cbrB_crc_df["diff z-score (cbrB-crc)"] = (
    cbrB_crc_df["Z score_cbrB"] - cbrB_crc_df["Z score_crc"]
)
cbrB_crc_df.head()

# Label genes in input list
cbrB_crc_df["gene group"] = "none"
cbrB_crc_df.loc[argR_genes, "gene group"] = "ArgR regulon"
cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "cbrB", "gene group"] = "cbrB"
cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "crc", "gene group"] = "crc"

print(cbrB_crc_df.shape)
cbrB_crc_df.head()

# +
# Drop cbrB and crc rows for plot
# Note: The knockout genes (cbrB and crc) were removed from the volcano plots below since
# they are "outliers" and the scaling for the z-score could be illustrated more clearly.

drop_idx = cbrB_crc_df[
    (cbrB_crc_df["gene group"] == "cbrB") | (cbrB_crc_df["gene group"] == "crc")
].index


cbrB_crc_df = cbrB_crc_df.drop(drop_idx)
cbrB_crc_df.shape

# +
# Plot
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["gene group"] == "none"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["lightgrey"],
    linewidth=0,
    legend=False,
    s=70,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["gene group"] == "ArgR regulon"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["red"],
    linewidth=0,
    s=70,
)

# Highlight specific genes mentioned in the paper (cbrB, crc, argA, aotJQMP)
"""h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "cbrB"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=1,
    palette=["#f2b4b4"],
    linewidth=0,
    legend=False,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "crc"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=1,
    palette=["#facf5a"],
    linewidth=0,
    legend=False,
)"""
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "argA"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["lightgrey"],
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "aotJ"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["red"],
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "aotQ"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["red"],
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "aotM"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["red"],
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
h = sns.scatterplot(
    data=cbrB_crc_df[cbrB_crc_df["Gene Name"] == "aotP"],
    x="logFC (Real)_cbrB",
    y="diff z-score (cbrB-crc)",
    hue="gene group",
    alpha=0.5,
    palette=["red"],
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)

# Add text labels for those genes mentioned above
"""plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "cbrB", "logFC (Real)_cbrB"] + 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "cbrB", "diff z-score (cbrB-crc)"]
    - 2,
    s="$cbrB$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "crc", "logFC (Real)_cbrB"] + 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "crc", "diff z-score (cbrB-crc)"] - 1,
    s="$crc$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "argA", "logFC (Real)_cbrB"] + 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "argA", "diff z-score (cbrB-crc)"],
    s="$argA$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotJ", "logFC (Real)_cbrB"] - 0.4,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotJ", "diff z-score (cbrB-crc)"]
    + 4,
    s="$aotJ$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotQ", "logFC (Real)_cbrB"] + 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotQ", "diff z-score (cbrB-crc)"]
    - 2,
    s="$aotQ$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotM", "logFC (Real)_cbrB"] - 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotM", "diff z-score (cbrB-crc)"]
    + 2,
    s="$aotM$",
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotP", "logFC (Real)_cbrB"] - 0.5,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotP", "diff z-score (cbrB-crc)"]
    - 6,
    s="$aotP$",
)
"""
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "argA", "logFC (Real)_cbrB"],
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "argA", "diff z-score (cbrB-crc)"]
    + 1,
    s="$argA$",
    size=12,
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotJ", "logFC (Real)_cbrB"] - 0.7,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotJ", "diff z-score (cbrB-crc)"],
    s="$aotJ$",
    size=12,
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotQ", "logFC (Real)_cbrB"] + 0.1,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotQ", "diff z-score (cbrB-crc)"],
    s="$aotQ$",
    size=12,
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotM", "logFC (Real)_cbrB"],
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotM", "diff z-score (cbrB-crc)"]
    - 2,
    s="$aotM$",
    size=12,
)
plt.text(
    x=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotP", "logFC (Real)_cbrB"] - 0.7,
    y=cbrB_crc_df.loc[cbrB_crc_df["Gene Name"] == "aotP", "diff z-score (cbrB-crc)"]
    - 3,
    s="$aotP$",
    size=12,
)


# Add traditional thresholds
h.axhline(0.0, c="black", lw=0.7, ls="--")
h.axvline(1, c="black", lw=0.7, ls="--")
h.axvline(-1, c="black", lw=0.7, ls="--")

plt.xlim(-4, 4)
plt.ylim(-20, 20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Make axis thicker
for _, s in h.spines.items():
    s.set_linewidth(1.5)

# Move location of legend
plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0, fontsize=12)

h.set_xlabel(r"Log$_{2}$ FC", fontsize=16, fontname="Verdana")
h.set_ylabel("cbrB z-score - crc z-score", fontsize=16, fontname="Verdana")
# -

# Save
h.figure.savefig(
    "cbrB_crc_zscore_compare.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# ## Volcano plot with genes colored by z-score

# +
# Read template DE stats
template_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_cbrB_v_WT.tsv"
)
template_summary_df = pd.read_csv(
    template_summary_filename, sep="\t", index_col=0, header=0
)

# Take -log10 of adjusted p-value
template_summary_df["padj_log10"] = -np.log10(template_summary_df["Adj P-value (Real)"])
# -

# Drop cbrB and crc rows for plotting -- see comment above for reason
print(template_summary_df.shape)
template_summary_df.head()

template_summary_df = template_summary_df.drop(drop_idx)
print(template_summary_df.shape)

# +
# Plot
f1 = sns.scatterplot(
    data=template_summary_df,
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    palette="rocket_r",
    linewidth=0,
    legend=False,
    s=70,
)

# Highlight specific common DEGs mentioned in the text
f1 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "pqsA"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f1 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "nosZ"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f1 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "pqsE"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f1 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "ccoP2"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)

# Add text labels for those genes mentioned above
"""plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsA", "logFC (Real)"
    ]
    - 0.9,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsA", "padj_log10"]
    - 0.2,
    s="$pqsA$",
)"""
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsA", "logFC (Real)"
    ]
    - 0.7,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsA", "padj_log10"]
    - 0.2,
    s="$pqsA$",
    size=12,
)
"""plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "nosZ", "logFC (Real)"
    ]
    + 0.2,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "nosZ", "padj_log10"],
    s="$nosZ$",
)"""
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "nosZ", "logFC (Real)"
    ]
    + 0.1,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "nosZ", "padj_log10"],
    s="$nosZ$",
    size=12,
)
"""plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsE", "logFC (Real)"
    ]
    + 0.2,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsE", "padj_log10"],
    s="$pqsE$",
)"""
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsE", "logFC (Real)"
    ]
    + 0.15,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsE", "padj_log10"],
    s="$pqsE$",
    size=12,
)
"""plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "ccoP2", "logFC (Real)"
    ]
    + 0.2,
    y=template_summary_df.loc[
        template_summary_df["Gene Name"] == "ccoP2", "padj_log10"
    ],
    s="$ccoP2$",
)"""
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "ccoP2", "logFC (Real)"
    ]
    + 0.15,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "ccoP2", "padj_log10"]
    - 0.2,
    s="$ccoP2$",
    size=12,
)

# Add traditional thresholds
f1.axhline(-np.log10(0.05), c="black", lw=0.7, ls="--")
f1.axvline(1, c="black", lw=0.7, ls="--")
f1.axvline(-1, c="black", lw=0.7, ls="--")

f1.set_ylabel(r"-Log$_{10}$ (FDR adjusted p-value)", fontsize=16, fontname="Verdana")
f1.set_xlabel(r"-Log$_{2}$ FC (cbrB/WT)", fontsize=16, fontname="Verdana")
# f1.set_title("$cbrB$ vs WT", fontsize=14, fontname="Verdana")

plt.xlim(-4, 4)
plt.ylim(0, 4)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Make axis thicker
for _, s in f1.spines.items():
    s.set_linewidth(1.5)

norm = plt.Normalize(
    template_summary_df["Z score"].min(), template_summary_df["Z score"].max()
)
sm = plt.cm.ScalarMappable(cmap="rocket_r", norm=norm)
sm.set_array([])

cb = f1.figure.colorbar(sm)

cb.ax.tick_params(labelsize=12)
cb.set_label(label="z-score", size=12, fontname="Verdana")
# -

# Save
f1.figure.savefig(
    "cbrB_volcano_zscore_highlight.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# +
# Read template DE stats
template_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_crc_v_WT.tsv"
)
template_summary_df = pd.read_csv(
    template_summary_filename, sep="\t", index_col=0, header=0
)

# Take -log10 of adjusted p-value
template_summary_df["padj_log10"] = -np.log10(template_summary_df["Adj P-value (Real)"])

print(template_summary_df.shape)
# -

template_summary_df = template_summary_df.drop(drop_idx)
print(template_summary_df.shape)

# +
# Plot
f2 = sns.scatterplot(
    data=template_summary_df,
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    palette="rocket_r",
    linewidth=0,
    legend=False,
    s=70,
)

# Highlight specific common DEGs mentioned in the text
f2 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "pqsA"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f2 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "nosZ"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f2 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "pqsE"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)
f2 = sns.scatterplot(
    data=template_summary_df[template_summary_df["Gene Name"] == "ccoP2"],
    x="logFC (Real)",
    y="padj_log10",
    hue="Z score",
    alpha=0.5,
    edgecolor="black",
    linewidth=1,
    legend=False,
    s=70,
)

# Add text labels for those genes mentioned above
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsA", "logFC (Real)"
    ]
    - 1.4,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsA", "padj_log10"]
    - 0.15,
    s="$pqsA$",
    size=12,
)
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "nosZ", "logFC (Real)"
    ]
    + 0.2,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "nosZ", "padj_log10"],
    s="$nosZ$",
    size=12,
)
"""plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsE", "logFC (Real)"
    ],
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsE", "padj_log10"]
    - 0.3,
    s="$pqsE$",
)"""
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "pqsE", "logFC (Real)"
    ]
    - 1.2,
    y=template_summary_df.loc[template_summary_df["Gene Name"] == "pqsE", "padj_log10"]
    - 0.4,
    s="$pqsE$",
    size=12,
)
plt.text(
    x=template_summary_df.loc[
        template_summary_df["Gene Name"] == "ccoP2", "logFC (Real)"
    ]
    + 0.2,
    y=template_summary_df.loc[
        template_summary_df["Gene Name"] == "ccoP2", "padj_log10"
    ],
    s="$ccoP2$",
    size=12,
)

# Add traditional thresholds
f2.axhline(-np.log10(0.05), c="black", lw=0.7, ls="--")
f2.axvline(1, c="black", lw=0.7, ls="--")
f2.axvline(-1, c="black", lw=0.7, ls="--")

f2.set_ylabel(r"-Log$_{10}$ (FDR adjusted p-value)", fontsize=16, fontname="Verdana")
f2.set_xlabel(r"-Log$_{2}$ FC (crc/WT)", fontsize=16, fontname="Verdana")
# f2.set_title("$crc$ vs WT", fontsize=14, fontname="Verdana")

plt.xlim(-4, 6)
plt.ylim(-0, 5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Make axis thicker
for _, s in f2.spines.items():
    s.set_linewidth(1.5)

norm = plt.Normalize(
    template_summary_df["Z score"].min(), template_summary_df["Z score"].max()
)
sm = plt.cm.ScalarMappable(cmap="rocket_r", norm=norm)
sm.set_array([])
cb = f2.figure.colorbar(sm)

cb.ax.tick_params(labelsize=12)
cb.set_label(label="z-score", size=12, fontname="Verdana")
# -

# Save
f2.figure.savefig(
    "crc_volcano_zscore_highlight.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
