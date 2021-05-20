"""
Author: Alexandra Lee
Date Created: 18 December 2020

This script provide supporting functions to run analysis notebooks.

This script includes functions to plot data
"""
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def make_volcano_template_highlight_DEG(
    template_DE_stats_filename,
    project_id,
    pval_name,
    logFC_name,
    output_figure_filename,
):
    """
    This function creates volcano plot of template experiment
    highlighting traditional DEGs

    Arguments
    ----------
    template_DE_stats_filename: str
        File containing DE statistics for template experiment
    project_id: str
        Experiment identifier
    pval_name: "padj" or "adj.P.Val"
    logFC_name: "logFC" or "log2FoldChange"
    output_figure_filename: str
        File to save figure to
    """

    # Read template DE stats
    template_DE_stats_df = pd.read_csv(
        template_DE_stats_filename, sep="\t", index_col=0, header=0
    )

    # Take -log10 of adjusted p-value
    template_DE_stats_df["padj_log10"] = -np.log10(template_DE_stats_df[pval_name])

    # Label DEGs by traditional criteria
    # log2FC > 1
    # padj < 0.05
    template_DE_stats_df["gene group"] = "none"
    template_DE_stats_df.loc[
        (abs(template_DE_stats_df[logFC_name]) > 1)
        & (template_DE_stats_df[pval_name] < 0.05),
        "gene group",
    ] = "DEG"

    # Plot
    colors = ["lightgrey", "#a1dab4ff"]

    f = sns.scatterplot(
        data=template_DE_stats_df,
        x=logFC_name,
        y="padj_log10",
        hue="gene group",
        hue_order=["none", "DEG"],
        style="gene group",
        markers={"none": ".", "DEG": "o"},
        palette=colors,
        linewidth=0,
        alpha=0.5,
    )
    handles, labels = f.get_legend_handles_labels()
    f.legend([handles[1]], [labels[1]], loc="upper right")

    f.set_xlabel(r"log$_2$ Fold Change", fontsize=14, fontname="Verdana")
    f.set_ylabel(r"-log$_{10}$ (FDR adjusted p-value)", fontsize=14, fontname="Verdana")
    f.set_title(f"Template experiment ({project_id})", fontsize=16, fontname="Verdana")

    f.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


def make_volcano_simulated_highlight_DEG(
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    ncols,
    nrows,
    fig_width,
    fig_height,
    output_figure_filename,
):
    """
    This function makes multiple volcano plots of example simulated experiments
    and highlights traditional DEGs

    Arguments
    ----------
    template_DE_stats_filename: str
        File containing DE statistics for template experiment
    project_id: str
        Experiment identifier
    pval_name: "padj" or "adj.P.Val"
    logFC_name: "logFC" or "log2FoldChange"
    num_simulated: int
        Number of simulated experiments
    ncols: int
        Number of columns in facet plot
    nrows: int
        Number of rows in facet plot
    fig_width: int
        Width of figure
    fig_height: ing
        Height of figure
    output_figure_filename: str
        File to save figure to
    """
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(fig_width, fig_height))

    axes = axes.ravel()

    for i in range(num_simulated):

        # Get filename
        simulated_DE_stats_filename = os.path.join(
            simulated_DE_stats_dir, f"DE_stats_simulated_data_{project_id}_{i}.txt",
        )

        # Read simulated DE stats
        simulated_DE_stats_df = pd.read_csv(
            simulated_DE_stats_filename, sep="\t", index_col=0, header=0
        )

        # Take -log10 of adjusted p-value
        simulated_DE_stats_df["padj_log10"] = -np.log10(
            simulated_DE_stats_df[pval_name]
        )

        # Label DEGs by traditional criteria
        # log2FC > 1
        # padj < 0.05
        simulated_DE_stats_df["gene group"] = "none"
        simulated_DE_stats_df.loc[
            (abs(simulated_DE_stats_df[logFC_name]) > 1)
            & (simulated_DE_stats_df[pval_name] < 0.05),
            "gene group",
        ] = "DEG"

        # Plot
        colors = ["lightgrey", "#a1dab4ff"]

        f = sns.scatterplot(
            data=simulated_DE_stats_df,
            x=logFC_name,
            y="padj_log10",
            hue="gene group",
            hue_order=["none", "DEG"],
            style="gene group",
            markers={"none": ".", "DEG": "o"},
            palette=colors,
            linewidth=0,
            alpha=0.5,
            legend=("full" if i == 0 else False),
            ax=axes[i],
        )

        axes[i].set_ylabel("")
        axes[i].set_xlabel("")
        if i == 0:
            handles, labels = f.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right")
            f.legend_.remove()

    fig.text(0.5, 0.0, r"log$_2$ Fold Change", ha="center", fontsize=14, fontname="Verdana")
    fig.text(
        0.08,
        0.5,
        r"-log$_{10}$ (FDR adjusted p-value)",
        va="center",
        rotation="vertical",
        fontsize=14,
        fontname="Verdana",
    )
    fig.suptitle(
        f"Example simulated experiments based on {project_id}",
        fontsize=16,
        fontname="Verdana",
    )

    fig.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


def get_generic_genes(gene_summary_filename, generic_threshold):
    """
    This function returns generic gene ids

    Arguments
    ----------
    gene_summary_filename: str
        File containing summary statistics for each gene
    generic_threshold: int
        Rank threshold to use to identify generic genes
    """
    # Read summary file
    summary_data = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)

    # Find all genes above generic_threshold
    generic_gene_ids = list(
        summary_data[summary_data["Rank (simulated)"] >= generic_threshold].index
    )

    return generic_gene_ids


def get_specific_genes(gene_summary_filename, num_top_genes):
    """
    This function returns specific gene ids

    Arguments
    ----------
    gene_summary_filename: str
        File containing summary statistics for each gene
    num_top_genes: int
        Number of genes with highest z-score
    """
    # Read summary file
    summary_data = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)

    # Find all genes above generic_threshold
    specific_gene_ids = list(summary_data.nlargest(num_top_genes, "Z score").index)
    print(len(specific_gene_ids))

    return specific_gene_ids


def make_volcano_template_highlight_generic_specific(
    gene_summary_filename,
    generic_threshold,
    num_specific_genes,
    template_DE_stats_filename,
    project_id,
    pval_name,
    logFC_name,
    output_figure_filename,
):
    """
    This function creates volcano plot of template experiment
    highlighting generic and specific genes.

    Arguments
    ----------
    gene_summary_filename: str
        File containing summary statistics for each gene
    generic_threshold: int
        Rank threshold to use to identify generic genes
    num_specific_genes: int
        Number of top Z-scoring genes. These will be
        considered specific genes
    template_DE_stats_filename: str
        File containing DE statistics for template experiment
    project_id: str
        Experiment identifier
    pval_name: "padj" or "adj.P.Val"
    logFC_name: "logFC" or "log2FoldChange"
    output_figure_filename: str
        File to save figure to
    """
    # Get generic gene ids
    generic_gene_ids = get_generic_genes(gene_summary_filename, generic_threshold)

    # Get specific gene ids
    specific_gene_ids = get_specific_genes(gene_summary_filename, num_specific_genes)

    # Read template DE stats
    template_DE_stats_df = pd.read_csv(
        template_DE_stats_filename, sep="\t", index_col=0, header=0
    )

    # Take -log10 of adjusted p-value
    template_DE_stats_df["padj_log10"] = -np.log10(template_DE_stats_df[pval_name])

    # Label generic genes
    template_DE_stats_df["gene group"] = "none"
    template_DE_stats_df.loc[generic_gene_ids, "gene group"] = "generic"
    template_DE_stats_df.loc[specific_gene_ids, "gene group"] = "specific"

    # Plot
    colors = ["lightgrey", "#2c7fb8", "red"]

    f = sns.scatterplot(
        data=template_DE_stats_df,
        x=logFC_name,
        y="padj_log10",
        hue="gene group",
        hue_order=["none", "generic", "specific"],
        style="gene group",
        markers={"none": ".", "generic": ".", "specific": "o"},
        palette=colors,
        linewidth=0,
        alpha=0.5,
    )

    f.set_xlabel(r"log$_2$ Fold Change", fontsize=14, fontname="Verdana")
    f.set_ylabel(r"-log$_{10}$ (FDR adjusted p-value)", fontsize=14, fontname="Verdana")
    f.set_title(f"Template experiment ({project_id})", fontsize=16, fontname="Verdana")

    f.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


def make_volcano_simulated_highlight_generic_specific(
    gene_summary_filename,
    generic_threshold,
    num_specific_genes,
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    ncols,
    nrows,
    fig_width,
    fig_height,
    output_figure_filename,
):

    """
    This function makes multiple volcano plots of example simulated experiments
    and highlights generic and specific genes on all volcano plots

    Arguments
    ----------
    gene_summary_filename: str
        File containing summary statistics for each gene
    generic_threshold: int
        Rank threshold to use to identify generic genes
    num_specific_genes: int
        Number of top Z-scoring genes. These will be
        considered specific genes
    template_DE_stats_filename: str
        File containing DE statistics for template experiment
    project_id: str
        Experiment identifier
    pval_name: "padj" or "adj.P.Val"
    logFC_name: "logFC" or "log2FoldChange"
    num_simulated: int
        Number of simulated experiments
    ncols: int
        Number of columns in facet plot
    nrows: int
        Number of rows in facet plot
    fig_width: int
        Width of figure
    fig_height: ing
        Height of figure
    output_figure_filename: str
        File to save figure to
    """
    # Get generic gene ids
    generic_gene_ids = get_generic_genes(gene_summary_filename, generic_threshold)

    # Get specific gene ids
    specific_gene_ids = get_specific_genes(gene_summary_filename, num_specific_genes)

    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(fig_width, fig_height))

    axes = axes.ravel()

    for i in range(num_simulated):

        # Get filename
        simulated_DE_stats_filename = os.path.join(
            simulated_DE_stats_dir, f"DE_stats_simulated_data_{project_id}_{i}.txt",
        )

        # Read simulated DE stats
        simulated_DE_stats_df = pd.read_csv(
            simulated_DE_stats_filename, sep="\t", index_col=0, header=0
        )

        # Take -log10 of adjusted p-value
        simulated_DE_stats_df["padj_log10"] = -np.log10(
            simulated_DE_stats_df[pval_name]
        )

        # Label generic genes
        simulated_DE_stats_df["gene group"] = "none"
        simulated_DE_stats_df.loc[generic_gene_ids, "gene group"] = "generic"
        simulated_DE_stats_df.loc[specific_gene_ids, "gene group"] = "specific"

        ## TO DO:
        # Add threshold for logFC and pvalue?

        # Plot
        colors = ["lightgrey", "#2c7fb8", "red"]

        f = sns.scatterplot(
            data=simulated_DE_stats_df,
            x=logFC_name,
            y="padj_log10",
            hue="gene group",
            hue_order=["none", "generic", "specific"],
            style="gene group",
            markers={"none": ".", "generic": ".", "specific": "o"},
            palette=colors,
            linewidth=0,
            alpha=0.5,
            legend=("full" if i == 0 else False),
            ax=axes[i],
        )

        axes[i].set_ylabel("")
        axes[i].set_xlabel("")
        if i == 0:
            handles, labels = f.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right")
            f.legend_.remove()

    fig.text(0.5, 0.0, r"log$_2$ Fold Change", ha="center", fontsize=14, fontname="Verdana")
    fig.text(
        0.08,
        0.5,
        r"-log$_{10}$ (FDR adjusted p-value)",
        va="center",
        rotation="vertical",
        fontsize=14,
        fontname="Verdana",
    )
    fig.suptitle(
        f"Example simulated experiments based on {project_id}",
        fontsize=16,
        fontname="Verdana",
    )

    fig.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )
