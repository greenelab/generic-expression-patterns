
# coding: utf-8

# # Visualize gene expression
# 
# This notebook visualizes the gene expression data for the template and simulated experiments in order to:
# 1. Validate that the structure of the gene expression data and simulated data are consistent
# 2. To visualize the signal that is in the experiments

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import os
import pandas as pd
import umap
import pickle
import glob
import seaborn as sns
from sklearn.decomposition import PCA
from keras.models import load_model
import plotnine as pn

from ponyo import utils
from generic_expression_patterns_modules import plot


# ## Load config parameters

# In[3]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)

params = utils.read_config(config_filename)


# In[4]:


# Load config params

local_dir = params['local_dir']
project_id = params['project_id']
num_simulated = params['num_simulated']

pval_name = "adj.P.Val"
logFC_name = "logFC"
run=0

# Settings for running visualization using pseudomonas config file
vae_model_dir = os.path.join(base_dir,"pseudomonas_analysis", "models", "NN_2500_30")
template_filename = params['processed_template_filename']
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = params['scaler_filename']


# ## Volcano plots

# In[5]:


# Check number of DEGs
template_DE_stats_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_template_data_{project_id}_real.txt"
)

template_DE_stats = pd.read_csv(
    template_DE_stats_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

selected = template_DE_stats[(template_DE_stats[pval_name]<0.01) & (abs(template_DE_stats[logFC_name])>1)]
print(selected.shape)


# In[6]:


plot.make_volcano_plot_template(
    template_DE_stats_filename,
    project_id,
    pval_name,
    logFC_name
)


# In[7]:


simulated_DE_stats_dir = os.path.join(local_dir, "DE_stats")

plot.make_volcano_plot_simulated(
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    5,
    5,
    20,
    15
)


# ## Explore flat volcano plots
# 
# Why are some volcano plots flat?

# In[8]:


# Lets look at the expression data for the volcano plot with DEGs
simulated_expression_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_11.txt"
)

simulated_expression = pd.read_csv(
    simulated_expression_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

simulated_expression


# In[9]:


sns.distplot(simulated_expression.loc["GSM822708_wtLB_A.CEL"])
sns.distplot(simulated_expression.loc["GSM822712_delta_cbrBLB_A.CEL"])


# In[10]:


# Lets look at the expression data for the flat volcano plot
simulated_expression_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_5.txt"
)

flat_simulated_expression = pd.read_csv(
    simulated_expression_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

flat_simulated_expression


# In[11]:


sns.distplot(flat_simulated_expression.loc["GSM822708_wtLB_A.CEL"])
sns.distplot(flat_simulated_expression.loc["GSM822712_delta_cbrBLB_A.CEL"])


# In[12]:


# Lets look at the DE stats associated with this flat volcano plot
simulated_DE_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_simulated_data_{project_id}_5.txt"
)

flat_simulated_DE = pd.read_csv(
    simulated_DE_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

flat_simulated_DE


# In[13]:


# Lets look at the DE stats associated with DEGs volcano plot
simulated_DE_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_simulated_data_{project_id}_11.txt"
)

simulated_DE = pd.read_csv(
    simulated_DE_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

simulated_DE


# Flat volcano plots are due to samples between the two groups are very similar, so the log2 fold change is very small 10^-1 for each gene and so the adjusted p-values are nearly the same for each gene. Another artifact of the VAE variance shrinkage.

# ## Explore leveling out of volcano plots
# 
# Why do volcano plots level out at at then extend outward? Is this an issue with the plotting?

# In[14]:


import matplotlib.pyplot as plt    
def make_volcano_plot_simulated_notransform(
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    ncols,
    nrows,
    fig_width,
    fig_height,
):
    """
	This function makes multiple volcano plots of example simulated experiments

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
        simulated_DE_stats_df["padj_log10"] = -(simulated_DE_stats_df[pval_name])

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
        colors = ["lightgrey", "#2c7fb8"]

        if i == 0:
            f = sns.scatterplot(
                data=simulated_DE_stats_df,
                x=logFC_name,
                y="padj_log10",
                hue="gene group",
                hue_order=["none", "DEG"],
                style="gene group",
                markers={"none": ".", "DEG": "o",},
                palette=colors,
                linewidth=0,
                alpha=0.5,
                legend="full",
                ax=axes[i],
            )

            axes[i].set_ylabel("")
            axes[i].set_xlabel("")
            handles, labels = f.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right")
            f.legend_.remove()

        else:
            f = sns.scatterplot(
                data=simulated_DE_stats_df,
                x=logFC_name,
                y="padj_log10",
                hue="gene group",
                hue_order=["none", "DEG"],
                style="gene group",
                markers={"none": ".", "DEG": "o",},
                palette=colors,
                linewidth=0,
                alpha=0.5,
                legend=False,
                ax=axes[i],
            )

            axes[i].set_ylabel("")
            axes[i].set_xlabel("")

    fig.text(0.5, 0.0, "log2 Fold Change", ha="center", fontsize=14, fontname="Verdana")
    fig.text(
        0.08,
        0.5,
        "-(FDR adjusted p-value)",
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


# In[15]:


simulated_DE_stats_dir = os.path.join(local_dir, "DE_stats")

make_volcano_plot_simulated_notransform(
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    5,
    5,
    20,
    15
)


# There is still a leveling off but it doesn't look as dramatic if we don't take the log10 but just negate the adjusted p-values.

# In[25]:


# Let's look at the distribution of adjusted p-value scores for those volcano
# plots that level off as increase logFC
simulated_DE_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_simulated_data_{project_id}_1.txt"
)

leveloff_simulated_DE = pd.read_csv(
    simulated_DE_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)


# In[26]:


# Compare the distribution of adjusted p-value scores
# for the volcano plot that doesn't level off as much
# as we increase logFC
simulated_DE_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_simulated_data_{project_id}_2.txt"
)

simulated_DE = pd.read_csv(
    simulated_DE_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)


# In[27]:


sns.distplot(leveloff_simulated_DE[pval_name])
sns.distplot(simulated_DE[pval_name])


# Looks like there is a peak of adjusted p-values at the minimum range of the distribution which is causing the leveling out (i.e. there are many genes with a similar low adjusted p-value). I am guessing this is also a result of the VAE shrinkage, where instead of genes having varying logFC, genes are compressed such that there are groups of genes with similar logFC and therefore similar adjusted p-values.
