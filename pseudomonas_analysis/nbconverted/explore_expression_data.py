#!/usr/bin/env python
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


# Heatmap visualization
sns.heatmap(simulated_expression.iloc[:,10:20])


# In[11]:


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


# In[12]:


sns.distplot(flat_simulated_expression.loc["GSM822708_wtLB_A.CEL"])
sns.distplot(flat_simulated_expression.loc["GSM822712_delta_cbrBLB_A.CEL"])


# In[13]:


# Heatmap visualization
sns.heatmap(flat_simulated_expression.iloc[:,10:20])


# In[14]:


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


# In[15]:


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

# In[16]:


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

        # Take negative of adjusted p-value
        simulated_DE_stats_df["padj_neg"] = -(simulated_DE_stats_df[pval_name])

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
                y="padj_neg",
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
                y="padj_neg",
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


# In[17]:


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

# In[18]:


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


# In[19]:


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


# In[20]:


sns.distplot(leveloff_simulated_DE[pval_name])
sns.distplot(simulated_DE[pval_name])


# In[21]:


# Get genes in peak of distribution (blue)
gene_ids = list(leveloff_simulated_DE[leveloff_simulated_DE[pval_name]<0.4].index)
print(leveloff_simulated_DE[leveloff_simulated_DE[pval_name]<0.4].head(20))

leveloff_simulated_expression_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_1.txt"
)

leveloff_simulated_expression = pd.read_csv(
    leveloff_simulated_expression_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)
sns.heatmap(leveloff_simulated_expression[gene_ids[0:20]])


# In[22]:


# Get genes in peak of distribution (orange)
gene_ids = list(simulated_DE[simulated_DE[pval_name]<0.4].index)
print(simulated_DE[simulated_DE[pval_name]<0.4].head(20))

simulated_expression_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_2.txt"
)

simulated_expression = pd.read_csv(
    simulated_expression_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)
sns.heatmap(simulated_expression[gene_ids[0:20]])


# In[23]:


# Let's look at the distribution of adjusted p-value scores for those volcano
# plots that is completely flat
simulated_DE_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_simulated_data_{project_id}_3.txt"
)

flat_simulated_DE = pd.read_csv(
    simulated_DE_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)


# In[24]:


sns.distplot(flat_simulated_DE[pval_name])


# In[25]:


sns.distplot(simulated_DE[pval_name])


# In[26]:


# Get genes in peak of distribution (flat)
gene_ids = list(flat_simulated_DE.index)
print(flat_simulated_DE.head(20))

flat_simulated_expression_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_3.txt"
)

flat_simulated_expression = pd.read_csv(
    flat_simulated_expression_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)
sns.heatmap(flat_simulated_expression[gene_ids[0:20]])


# Looks like there is a peak of adjusted p-values at the minimum range of the distribution which is causing the leveling out (i.e. there are many genes with a similar low adjusted p-value). I am guessing this is also a result of the VAE shrinkage, where instead of genes having varying logFC, genes are compressed such that there are groups of genes with similar logFC and therefore similar adjusted p-values.
# 
# The first heatmap shows the expression in the simulated experiment with flat top (blue) for those genes with low adjusted p-values (i.e. in the peak of the distribution). The second heatmap shows the expression of the simulated experiment with a more V-shape (orange) for those genes with low adjusted p-values (i.e. those with adjusted p-values = 0-0.4). For the more V-shaped experiment, it looks there was more consistency amongst samples within the group (i.e. WT vs mutant). Depending on the simulation experiment generated, there might be some noise created by the VAE.
# 
# A similar trend is seen for the volcano plot that is completely flat as well. There is more noise in the gene expression data.

# ## Expression of template experiment
# 
# Do the samples have a clear separation in gene space (WT vs mutant)?

# In[27]:


normalized_compendium_data = pd.read_csv(normalized_compendium_filename, sep="\t", index_col=0, header=0)
template_data = pd.read_csv(template_filename, sep="\t", index_col=0, header=0)


# In[28]:


print(template_data.shape)
template_data


# In[29]:


# If template experiment included in training compendium
# Get normalized template data
sample_ids = list(template_data.index)
normalized_template_data = normalized_compendium_data.loc[sample_ids]

print(normalized_template_data.shape)
normalized_template_data.head()


# In[30]:


# Label samples 
wt_sample_ids = ["GSM822708_wtLB_A.CEL", "GSM822709_wtLB_B.CEL"]
mutant_sample_ids = ["GSM822712_delta_cbrBLB_A.CEL", "GSM822713_delta_cbrBLB_B.CEL"]
normalized_compendium_data['sample group'] = "compendium"
normalized_template_data.loc[wt_sample_ids, 'sample group'] = "template_WT"
normalized_template_data.loc[mutant_sample_ids, 'sample group'] = "template_mutant"


# In[31]:


normalized_all_data = pd.concat([normalized_template_data,
                                 normalized_compendium_data
])


# In[32]:


# Plot

# Drop label column
normalized_all_data_numeric = normalized_all_data.drop(['sample group'], axis=1)

model = umap.UMAP(random_state=1).fit(normalized_all_data_numeric)

normalized_all_data_UMAPencoded = model.transform(normalized_all_data_numeric)
normalized_all_data_UMAPencoded_df = pd.DataFrame(data=normalized_all_data_UMAPencoded,
                                         index=normalized_all_data.index,
                                         columns=['1','2'])

# Add back label column
normalized_all_data_UMAPencoded_df['sample group'] = normalized_all_data['sample group']

# Plot
fig = pn.ggplot(normalized_all_data_UMAPencoded_df, pn.aes(x='1', y='2'))
fig += pn.geom_point(pn.aes(color='sample group'), alpha=0.4)
fig += pn.labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'Gene expression data in gene space')
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align = "center",
    plot_background=pn.element_rect(fill='white'),
    legend_key=pn.element_rect(fill='white', colour='white'), 
    legend_title=pn.element_text(family='sans-serif', size=15),
    legend_text=pn.element_text(family='sans-serif', size=12),
    plot_title=pn.element_text(family='sans-serif', size=15),
    axis_text=pn.element_text(family='sans-serif', size=12),
    axis_title=pn.element_text(family='sans-serif', size=15)
    )
fig += pn.scale_color_manual(['#bdbdbd', 'red', 'blue'])
fig += pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1}))

fig += pn.scales.xlim(9,10)
print(fig)


# Based on a UMAP of the normalized gene expression data, it looks like there isn't a clear separation between WT and mutant samples, though there are only 2 samples per group so this type of clustering observation is limited.
# 
# **Takeaway:**
# 
# In trying to understand why there are these flat-tops to some of the volcano plots and why some volcano plots are completely flat, we found:
# 1. This behavior is _not_ a result of how we are plotting in python (there was some speculation about there being an issue with the numpy library used)
# 2. The latent space shifting we're doing seems to roughly preserve differences between groups (as seen in [this notebook](https://github.com/greenelab/simulate-expression-compendia/blob/master/Pseudo_experiments/create_heatmap.ipynb) where the structure of the samples is preserved but there is a different set of related genes that are DE. More information can be found in Figure 3D in [this paper](https://academic.oup.com/gigascience/article/9/11/giaa117/5952607)), but this signal can be muddled/noisy depending on where the experiment was shifted to (i.e. the representation that is found in that location can cause the experiment to have a more compressed difference between groups) as seen in the heatmaps. The heatmap of the two simulation experiments shows that some experiments have a more noisey distinction between groups (WT vs mutant) whereas the other simulation experiment has a more distinct difference where the within grouping is cleaner. This definitely points to the need to understand how this simulation process is working and how biology is represented in the latent space. This will definitely be a project for the future. For now we at least have an explanation for why we are observing these shapes in the volcano plots
