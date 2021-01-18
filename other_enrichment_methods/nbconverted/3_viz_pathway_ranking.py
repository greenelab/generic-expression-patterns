
# coding: utf-8

# # Visualize pathway ranking
# 
# This notebook will visualize pathway ranking obtained by the different enrichment analysis methods.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import numpy as np
import plotnine as pn

from ponyo import utils

np.random.seed(123)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)


# In[3]:


# Load params
dataset_name = params["dataset_name"]
project_id = params["project_id"]


# In[4]:


# Create dictionary of enrichment method: statistic
method_stats_dict = {
    "GSEA": "padj",
    "GSVA": "ES",
    "ROAST": "FDR",
    "CAMERA": "FDR",
    "ORA": "p.adjust"
}


# ## Get pathway summary data

# In[5]:


# Pathway summary files
gsea_pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}.tsv"
)
gsva_pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}_GSVA.tsv"
)
roast_pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}_ROAST.tsv"
)
camera_pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}_CAMERA.tsv"
)
ora_pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}_ORA.tsv"
)


# In[6]:


# Pathway summary data
gsea_pathway_summary = pd.read_csv(gsea_pathway_summary_filename, sep="\t", index_col=0, header=0)
gsva_pathway_summary = pd.read_csv(gsva_pathway_summary_filename, sep="\t", index_col=0, header=0)
roast_pathway_summary = pd.read_csv(roast_pathway_summary_filename, sep="\t", index_col=0, header=0)
camera_pathway_summary = pd.read_csv(camera_pathway_summary_filename, sep="\t", index_col=0, header=0)
ora_pathway_summary = pd.read_csv(ora_pathway_summary_filename, sep="\t", index_col=0, header=0)


# ## Format data for plotting

# In[7]:


gsea_pathway_summary.head()


# In[8]:


gsva_pathway_summary.head()


# In[9]:


roast_pathway_summary.head()


# In[10]:


camera_pathway_summary.head()


# In[11]:


ora_pathway_summary.head()


# ## Plot

# In[12]:


# define plotting function
def plot_significance_vs_ranking(summary_df, method_name, x_label, output_figure_filename):
    # Format input dataframe
    plot_df = pd.DataFrame(
        data={"Test statistic": summary_df[method_stats_dict[method_name]+" (Real)"].values,
          "Percentile rank": summary_df["Rank (simulated)"].rank(pct=True).values
         },
        index=summary_df.index
    )
    
    fig = pn.ggplot(plot_df, pn.aes(x='Test statistic', y='Percentile rank'))
    fig += pn.geom_point()
    fig += pn.geom_point(plot_df[plot_df['Percentile rank']>0.9],
                         pn.aes(x='Test statistic', y='Percentile rank'),
                         color='red'
                        )
    fig += pn.geom_text(pn.aes(label=
                               [x if plot_df.loc[x,'Percentile rank']>0.9 else "" for x in plot_df.index]),
                        ha='left',
                        va='top',
                        size=5
                       )
    fig += pn.labs(x = x_label,
                y = 'Percentile of ranking',
                title = f'{method_name} pathway statistics vs ranking')
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

    print(fig)
    
    # Save figure
    fig.save(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


# In[13]:


plot_significance_vs_ranking(gsea_pathway_summary, "GSEA", "adjusted p-value (BH)", "GSEA_pathway_ranking.svg")


# In[14]:


plot_significance_vs_ranking(gsva_pathway_summary, "GSVA", "Enrichment score", "GSVA_pathway_ranking.svg")


# In[15]:


plot_significance_vs_ranking(roast_pathway_summary, "ROAST", "FDR (BH)", "ROAST_pathway_ranking.svg")


# In[16]:


plot_significance_vs_ranking(camera_pathway_summary, "CAMERA", "FDR (BH)", "CAMERA_pathway_ranking.svg")


# In[17]:


plot_significance_vs_ranking(ora_pathway_summary, "ORA", "adjusted p-value (BH)", "ORA_pathway_ranking.svg")


# **Takeaway:**
# * Here are the results demonstrating that different enrichment methods can easily be plugged into our simulation workflow to identify generic gene sets
# * Overall, it appears that the most enrichment gene sets are also those that are most commonly found to be enriched. However, depending on the enrichment method, generic gene sets will vary slightly due to the different assumptions and modeling procedures. More details about the methods can be found in the [previous notebook](2_apply_enrichment_method.ipynb)
# * ROAST performs a self-contained test, which assesses the relevance of an individual biological process to the experiment at hand without reference to other genes in the genome. This lack of reference might account for why the scores are uniform across the different gene sets.
