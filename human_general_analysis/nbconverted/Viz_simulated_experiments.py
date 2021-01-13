
# coding: utf-8

# ## Visualize simulated experiment
# 
# This notebook will make volcano plots of 3 representative simulated experiments to demonstrate what the simulation approach creates from a biological perspective (i.e. experiments created contain a similar experiment structure but different biological patterns)

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from ponyo import utils


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)


# In[3]:


# Load params
local_dir = params["local_dir"]
project_id = params['project_id']

simulated_DE_stats_dir = os.path.join(
    local_dir, 
    "DE_stats"
)
template_DE_stats_filename = os.path.join(
    local_dir, 
    "DE_stats",
    f"DE_stats_template_data_{project_id}_real.txt"
    
)


# In[4]:


# Plot volcano plots
def make_volcano_plot_template(
    template_DE_stats_filename,
    project_id,
    out_filename
):
    
    # Read template DE stats
    template_DE_stats_df = pd.read_csv(
        template_DE_stats_filename,
        sep="\t",
        index_col=0,
        header=0
    )
    
    # Take -log10 of adjusted p-value
    template_DE_stats_df["padj_log10"] = -np.log10(template_DE_stats_df["padj"])

    # Label DEGs by traditional criteria
    # log2FC > 1
    # padj < 0.05
    template_DE_stats_df["gene group"] = "none"
    template_DE_stats_df.loc[(abs(template_DE_stats_df["log2FoldChange"])>1) &
                          (template_DE_stats_df["padj"] <0.05),
                              "gene group"
                         ] = "DEG"

    # Plot
    colors = ["lightgrey", "#2c7fb8"]

    f = sns.scatterplot(
       data=template_DE_stats_df,
        x="log2FoldChange",
        y="padj_log10",
        hue="gene group",
        hue_order=["none", "DEG"],
        style="gene group",
        markers={
            "none": ".",
            "DEG": "o",
        },
        palette=colors,
        linewidth=0,
        alpha=0.5,
        )
    
    f.set_xlabel("log2 Fold Change", fontsize=14, fontname="Verdana")
    f.set_ylabel("-log10(FDR adjusted p-value)", fontsize=14, fontname="Verdana")
    f.set_title(f"Template experiment ({project_id})", fontsize=16, fontname="Verdana")
    
    # Save plot
    f.figure.savefig(
        out_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )
        
    
def make_volcano_plot_simulated(
    simulated_DE_stats_dir,
    project_id,
    num_examples,
    out_filename
):
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))
    
    for i in range(num_examples):
        
        # Get filename
        simulated_DE_stats_filename = os.path.join(
            simulated_DE_stats_dir,
            f"DE_stats_simulated_data_{project_id}_{i}.txt"
        )
        
        # Read simulated DE stats
        simulated_DE_stats_df = pd.read_csv(
            simulated_DE_stats_filename,
            sep="\t",
            index_col=0,
            header=0
        )

        # Take -log10 of adjusted p-value
        simulated_DE_stats_df["padj_log10"] = -np.log10(simulated_DE_stats_df["padj"])

        # Label DEGs by traditional criteria
        # log2FC > 1
        # padj < 0.05
        simulated_DE_stats_df["gene group"] = "other gene"
        simulated_DE_stats_df.loc[(abs(simulated_DE_stats_df["log2FoldChange"])>1) &
                              (simulated_DE_stats_df["padj"] <0.05),
                                  "gene group"
                             ] = "DEG"
        
        # Plot
        colors = ["lightgrey", "#2c7fb8"]
        
        if i == 0:
            f = sns.scatterplot(
                data=simulated_DE_stats_df,
                x="log2FoldChange",
                y="padj_log10",
                hue="gene group",
                hue_order=["other gene", "DEG"],
                style="gene group",
                markers={"other gene": ".", "DEG": "o",},
                palette=colors,
                linewidth=0,
                alpha=0.5,
                legend="full",
                ax=axes[i],
            )

            axes[i].set_ylabel("")
            axes[i].set_xlabel("")
            
            # Note: We are creating a single global legend that apply
            # to all the facets of this figure. To do this using
            # matplotlib, we need to be a little creative here
            # and add the legend to a new location that is applied
            # to the figure and then remove the legend from the facet.
            handles, labels = f.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right")
            f.legend_.remove()

        else:
            f = sns.scatterplot(
                data=simulated_DE_stats_df,
                x="log2FoldChange",
                y="padj_log10",
                hue="gene group",
                hue_order=["other gene", "DEG"],
                style="gene group",
                markers={"other gene": ".", "DEG": "o",},
                palette=colors,
                linewidth=0,
                alpha=0.5,
                legend=False,
                ax=axes[i],
            )

            axes[i].set_ylabel("")
            axes[i].set_xlabel("")

    fig.text(0.5, 0.0, "log2 Fold Change",ha="center", fontsize=14, fontname="Verdana")
    fig.text(0.08, 0.5, "-log10(FDR adjusted p-value)", va="center", rotation="vertical", fontsize=14, fontname="Verdana")
    fig.suptitle(f"Example simulated experiments based on {project_id}", fontsize=16, fontname="Verdana")
    
    # Save plot
    fig.savefig(
        out_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


# In[5]:


# Plot template volcano
out_filename = "template_volcano.svg"
make_volcano_plot_template(template_DE_stats_filename, project_id, out_filename)


# In[6]:


# Plot simulated volcanos
out_filename = "example_simulated_volcano.svg"
make_volcano_plot_simulated(simulated_DE_stats_dir, project_id, 3, out_filename)


# **Takeaway:**
# * Simulated experiments have DEGs with a similar distribution of effect size changes (i.e. our simulation preserve the experimental design, so the "changes" are maintained)
# * However, most DEGs in the simulated experiments are down-regulated compared to the template experiment where genes are equally split between up-regulated and down-regulated. This indicates the different biological patterns that are generated by the simulated (i.e. the simulated experiments are different experiments compared to the template). Furthermore, there are inter-simulated experiments differences in the volcano plot, indicating that these are distinct new experiments. 
# * In both template and simulated experiments, DEGs are determines by p-value more than the log2 fold change (i.e minimally changed genes are signficant). This might suggest that very few genes were changed in this experiment? (Need to think about this more)
