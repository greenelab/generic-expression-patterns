
# coding: utf-8

# # Prepare reference gene file
# 
# This notebook creates a reference gene ranking file to use to compare SOPHIE generated gene ranking. The reference ranking information is obtained from [this repository](https://github.com/DartmouthStantonLab/GAPE). This [RDS object](https://github.com/DartmouthStantonLab/GAPE/blob/main/Pa_GPL84_refine_ANOVA_List_unzip.rds) contains 73 experiments. For each experiment, we will identify DEGs using log2FC > 1 and FDR < 0.05. We will rank genes by the proportion that they appeared as DE.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ponyo import utils


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

config_filename = os.path.abspath(os.path.join(base_dir,
                                           "configs",
                                           "config_pseudomonas_33245.tsv"))
params = utils.read_config(config_filename)


# In[3]:


# Load params
local_dir = params["local_dir"]
reference_gene_filename = os.path.join(local_dir,"Pa_GPL84_refine_ANOVA_List_unzip.rds")


# In[4]:


readRDS = ro.r['readRDS']


# In[5]:


reference_Robject = readRDS(reference_gene_filename)


# In[6]:


# For each experiment get df
# For each df, if label gene as DEGs based on log2FC>1 and FDR<0.05
# Concatenate series
num_experiments = len(reference_Robject)
reference_stats_df = pd.DataFrame()

for i in range(num_experiments):
    # Get df for experiment
    reference_df = pandas2ri.ri2py(reference_Robject[i])
    
    reference_df = reference_df.set_index("ProbeID")
    
    # Find DEGs
    degs_ids = list(reference_df[(abs(reference_df["Log2FC"])>1) & (reference_df["FDR"]<0.05)].index)
    reference_df["DEG"] = 0
    reference_df.loc[degs_ids, "DEG"] = 1
    
    # Create df with labels for if gene is DE or not
    if i == 0:
        reference_stats_df = reference_df["DEG"].to_frame("experiment_0") 
    else:
        reference_stats_df = pd.merge(
            reference_stats_df,
            reference_df["DEG"].to_frame(f"experiment_{i}"),
            left_index=True, 
            right_index=True,
            how='left'
        )
reference_stats_df 


# In[7]:


# Map `ProbeID` to `IntergenicSpotID` that contains PA#### IDs
example_reference_df = pandas2ri.ri2py(reference_Robject[0])
example_reference_df.set_index("ProbeID", inplace=True)

merged_df = pd.merge(reference_stats_df, example_reference_df, left_index=True, right_index=True)

# Get relevant columns (`IntergenicSpotID` and `experiment_*`)
experiment_cols = [col for col in merged_df.columns if 'experiment_' in col]
merged_df = merged_df[experiment_cols]

# Aggregate to get ranking of genes
merged_df["prop DEGs"] = merged_df.sum(axis=1)/num_experiments


# In[8]:


# Extract PA#### ids from `ProbeID`
# This will be used to compare against SOPHIE ranked genes
pao1_ids = [str_ls[0] for str_ls in merged_df.index.str.split("_")]
merged_df["gene id"] = pao1_ids


# In[9]:


merged_df.head()


# In[10]:


# Save file
# Here are the names that we will use for the comparison in notebook 2_identify_generic_genes_pathways.ipynb
# DE_prior_filename = output_filename
# ref_gene_col = "gene id"
# ref_rank_col = "prop DEGs"
merged_df.to_csv('GAPE_proportions.txt', sep="\t")

