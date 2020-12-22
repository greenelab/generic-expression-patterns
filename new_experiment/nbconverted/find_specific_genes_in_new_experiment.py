
# coding: utf-8

# # Application: new experiment
# 
# This notebook allows users to find specific genes in their experiment of interest using an existing VAE model
# 
# This notebook will generate a `generic_gene_summary_<experiment id>.tsv` file that contains a z-score per gene that indicates how specific a gene is the experiment in question.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from ponyo import utils
from generic_expression_patterns_modules import process, new_experiment_process, stats, ranking


# ## User input
# 
# User needs to define the following in the [config file](../configs/config_new_experiment.tsv):
# 
# 1. Template experiment. This is the experiment you are interested in studying
# 2. Training compendium used to train VAE, including unnormalized gene mapped version and normalized version
# 3. Scaler transform used to normalize the training compendium
# 4. Directory containing trained VAE model
# 5. Experiment id to label newly create simulated experiments
# 
# The user also needs to provide metadata files:
# 1. `<experiment id>_process_samples.tsv` contains 2 columns (sample ids, label that indicates if the sample is kept or removed). See [example](data/metadata/cis-gem-par-KU1919_process_samples.tsv). **Note: This file is not required if the user wishes to use all the samples in the template experiment file.**
# 2. `<experiment id>_groups.tsv` contains 2 columns: sample ids, group label to perform DE analysis. See [example](data/metadata/cis-gem-par-KU1919_groups.tsv)

# In[3]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_new_experiment.tsv")
)

params = utils.read_config(config_filename)


# In[4]:


# Load config params

# Local directory to store intermediate files
local_dir = params['local_dir']

# Number of simulated experiments to generate
num_runs = params['num_simulated']

# Directory containing trained VAE model
vae_model_dir = params['vae_model_dir']

# Dimension of latent space used in VAE model
latent_dim = params['latent_dim']

# ID for template experiment
# This ID will be used to label new simulated experiments
project_id = params['project_id']

# Template experiment filename
template_filename = params['raw_template_filename']
processed_template_filename = params['processed_template_filename']

# Training dataset used for existing VAE model
mapped_compendium_filename = params['mapped_compendium_filename']

# Normalized compendium filename
normalized_compendium_filename = params['normalized_compendium_filename']

# Scaler transform used to scale compendium data into 0-1 range for training
scaler_filename = params['scaler_filename']

# Test statistic used to rank genes by
col_to_rank_genes = params['rank_genes_by']


# In[5]:


# Load metadata files

# Load metadata file with processing information
sample_id_metadata_filename = os.path.join(
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    "data",
    "metadata",
    f"{project_id}_groups.tsv"
)


# In[6]:


# Output filename
gene_summary_filename = f"generic_gene_summary_{project_id}.tsv"


# ## Map template experiment to same space as compendium
# 
# TO DO: Consolidate these functions

# In[7]:


# Template experiment needs to be of the form sample x gene
transposed_template_filename = "/home/alexandra/Documents/Data/Generic_expression_patterns/Costello_BladderCancer_ResistantCells_Counts_12-8-20_transposed.txt"

new_experiment_process.transpose_save(template_filename, transposed_template_filename)


# In[8]:


# Check that the feature space matches between template experiment and VAE model.  
# (i.e. ensure genes in template and VAE model are the same).
mapped_template_experiment = new_experiment_process.compare_match_features(
    transposed_template_filename,
    mapped_compendium_filename
)
mapped_template_filename = transposed_template_filename


# In[9]:


# Scale template experiment to be within the same range as the
# normalized training dataset used for the VAE model
new_experiment_process.normalize_template_experiment(
    mapped_template_experiment,
    scaler_filename,
    processed_template_filename
)


# ## Simulate experiments based on template experiment
# 
# Embed template experiment into learned latent space and linearly shift template experiment to different locations of the latent space to create new experiments

# In[10]:


# Simulate experiments based on template experiment
normalized_data = pd.read_csv(normalized_compendium_filename, sep="\t", index_col=0, header=0)
processed_template_data = pd.read_csv(processed_template_filename, sep="\t", index_col=0, header=0)

for run_id in range(num_runs):
    new_experiment_process.embed_shift_template_experiment(
        normalized_data,
        processed_template_data,
        vae_model_dir,
        project_id,
        scaler_filename,
        local_dir,
        latent_dim,
        run_id
    )


# ## Process template and simulated experiments
# 
# * Remove samples not required for comparison
# * Make sure ordering of samples matches metadata for proper comparison
# * Make sure values are cast as integers if using DESeq
# * Filter lowly expressed genes if using DESeq

# In[11]:


if os.path.exists(sample_id_metadata_filename):
    stats.process_samples_for_DESeq(
        mapped_template_filename,
        sample_id_metadata_filename,
        metadata_filename,
        None,
        processed_template_filename
    )
    
    for i in range(num_runs):
        simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt"
        )
        stats.process_samples_for_DESeq(
        simulated_filename,
        sample_id_metadata_filename,
        metadata_filename
    )


# ## Differential expression analysis
# 
# TO DO: Add conditonal for which method limma vs DESeq
# 
# * If data is RNA-seq then use DESeq2 (using human_general_analysis model)
# * If data is microarray then use Limma (using human_cancer_analysis, pseudomonas_analysis models)

# In[12]:


# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)


# In[13]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/DE_analysis.R\'))\n\n# File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"\nget_DE_stats_DESeq(metadata_filename,\n                   project_id, \n                   processed_template_filename,\n                   "template",\n                   local_dir,\n                   "real")')


# In[14]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/DE_analysis.R\'))\n\n# Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_<project_id>_<n>.txt"\nfor (i in 0:(num_runs-1)){\n    simulated_data_filename <- paste(local_dir, \n                                     "pseudo_experiment/selected_simulated_data_",\n                                     project_id,\n                                     "_", \n                                     i,\n                                     ".txt",\n                                     sep = "")\n    \n    get_DE_stats_DESeq(metadata_filename,\n                       project_id, \n                       simulated_data_filename,\n                       "simulated",\n                       local_dir,\n                       i)\n}')


# ## Quick validation
# 
# Examine volcano plot of the template experiment and example simulated experiments.
# 
# TO DO: move this to new notebook

# In[15]:


# Quick validation
def make_volcano_plot_template(
    template_DE_stats_filename,
    project_id,
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
    
def make_volcano_plot_simulated(
    simulated_DE_stats_dir,
    project_id,
    num_examples
):
    fig, axes = plt.subplots(ncols=num_examples, nrows=1, figsize=(15, 4))
    
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
        simulated_DE_stats_df["gene group"] = "none"
        simulated_DE_stats_df.loc[(abs(simulated_DE_stats_df["log2FoldChange"])>1) &
                              (simulated_DE_stats_df["padj"] <0.05),
                                  "gene group"
                             ] = "DEG"
        
        # Plot
        colors = ["lightgrey", "#2c7fb8"]
        
        f = sns.scatterplot(
           data=simulated_DE_stats_df,
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
            legend=False,
            ax=axes[i],
            )
        
        axes[i].set_ylabel("")
        axes[i].set_xlabel("")
        
        
    fig.legend(labels=["DEGs", "other genes"], loc='center right')
    fig.text(0.5, 0.0, "log2 Fold Change",ha="center", fontsize=14, fontname="Verdana")
    fig.text(0.08, 0.5, "-log10(FDR adjusted p-value)", va="center", rotation="vertical", fontsize=14, fontname="Verdana")
    fig.suptitle(f"Example simulated experiments based on {project_id}", fontsize=16, fontname="Verdana")


# In[16]:


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

selected = template_DE_stats[(template_DE_stats['padj']<0.01) & (abs(template_DE_stats['log2FoldChange'])>1)]
print(selected.shape)


# In[17]:


make_volcano_plot_template(
    template_DE_stats_filename,
    project_id,
)


# In[18]:


simulated_DE_stats_dir = os.path.join(local_dir, "DE_stats")
num_examples = 3
make_volcano_plot_simulated(
    simulated_DE_stats_dir,
    project_id,
    num_examples
)


# ## Rank genes
# 
# Add description

# In[19]:


analysis_type = "DE"
template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
    template_DE_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_genes,
)


# ## Summary table
# 
# TO DO: Description of table columns
# 
# Note: If using DESeq, genes with NaN in `Adj P-value (Real)` column are those genes flagged because of the `cooksCutoff` parameter. The cook's distance as a diagnostic to tell if a single sample has a count which has a disproportionate impact on the log fold change and p-values. These genes are flagged with an NA in the pvalue and padj columns of the result table. For more information you can read [DESeq FAQs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)

# In[20]:


# Get summary table
summary_gene_ranks = ranking.generate_summary_table(
    template_DE_stats_filename,
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir,
    'gene',
    params
)

summary_gene_ranks.sort_values(by="Z score", ascending=False).head(10)


# In[21]:


summary_gene_ranks.isna().any()


# In[22]:


summary_gene_ranks[summary_gene_ranks.isna().any(axis=1)]


# In[23]:


# Save
summary_gene_ranks.to_csv(gene_summary_filename, sep='\t')

