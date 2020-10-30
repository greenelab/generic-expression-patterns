#!/usr/bin/env python
# coding: utf-8

# # Multiplier analysis
# 
# The goal of this notebook is to examine why genes were found to be generic. Specifically, this notebook is trying to answer the question: Are generic genes found in more multiplier latent variables compared to specific genes?
# 
# The PLIER model performs a matrix factorization of gene expression data to get two matrices: loadings (Z) and latent matrix (B). The loadings (Z) are constrained to aligned with curated pathways and gene sets specified by prior knowledge [Figure 1B of Traoni et. al.](). This ensure that some but not all latent variables capture known biology. The way PLIER does this is by applying a penalty such that the individual latent variables represent a few gene sets in order to make the latent variables more interpretable. Ideally there would be one latent variable associated with one gene set unambiguously.
# 
# While the PLIER model was trained on specific datasets, MULTIPLIER extended this approach to all of recount2, where the latent variables should correspond to specific pathways or gene sets of interest. Therefore, we will look at the coverage of generic genes versus specific genes across these MULTIPLIER latent variables.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import re
import pandas as pd

from generic_expression_patterns_modules import process


# In[2]:


# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "human_general_analysis")


# In[ ]:


# TO DO: Move all functions to file


# In[3]:


# Read in all files of the form "generic_gene_summary_*"
# For each file, return list of generic genes and list of specific genes
# Should it return dictionary of lists?
def get_gene_summary_files(data_dir):
    files = [os.path.join(data_dir,f) for f in os.listdir(data_dir) if re.match(r'generic_gene_summary_*', f)]
    return files

ls_data_files = get_gene_summary_files(data_dir)


# In[4]:


ls_data_files


# In[5]:


def get_generic_specific_genes(list_files, z_threshold):
    ls_genes = []
    for file in list_files:
        print(f"Reading data for {file}")
        data = pd.read_csv(file, sep="\t", index_col=0, header=0)
        print(data.shape)
        
        # Get predicted specific DEGs using z-score cutoff
        ls_specific_genes = list(
            (
                data[(data[f"Test statistic (Real)"] > 1)
                    & (data[f"abs(Z score)"] > z_threshold
                    )
                ]
                .set_index("Gene ID")
                .index
            )
        )
        print(f"No. of specific DEGs using z-score: {len(ls_specific_genes)}")

        # Get predicted generic DEGs using z-score cutoff
        ls_generic_genes = list(
            (
                data[
                    (data[f"Test statistic (Real)"] > 1)
                    & (data[f"abs(Z score)"]< z_threshold
                    )
                ]
                .set_index("Gene ID")
                .index
            )
        )
        print(f"No. of generic DEGs using z-score: {len(ls_generic_genes)}")
    
        ls_genes.append([ls_generic_genes, ls_specific_genes])
    
    return ls_genes

# TO DO: add more accurate description here
# Get predicted generic DEGs using z-score cutoff
# Z-score cutoff was found by calculating invnorm(0.05/17754). 
# To do this in python you can use the following code:
# from scipy.stats import norm
# norm.ppf((0.05/17754)/2)
# Here we are using a p-value = 0.05
# with a Bonferroni correction for 17754 tests, which are
# the number of P. aeruginosa genes

zscore_threshold = 4.68
ls_genes_out = get_generic_specific_genes(ls_data_files, zscore_threshold)

# TO DO:
# Is this how we want to define the genes? ranking?


# In[6]:


# Load multiplier models
# Converted formatted pickle files (loaded using phenoplier environment) from
# https://github.com/greenelab/phenoplier/blob/master/nbs/01_preprocessing/005-multiplier_recount2_models.ipynb
# into .tsv files
# Raw data was downloaded from https://figshare.com/articles/recount_rpkm_RData/5716033/4
multiplier_model_u = pd.read_csv("multiplier_model_u.tsv", sep="\t", index_col=0, header=0)
multiplier_model_z = pd.read_csv("multiplier_model_z.tsv", sep="\t", index_col=0, header=0)

print(multiplier_model_u.shape)
multiplier_model_u.head()


# In[7]:


print(multiplier_model_z.shape)
multiplier_model_z.head()


# In[17]:


# Get a rough sense for how many genes contribute to a given LV
# (i.e. how many genes have a value > 0 per LV)
(multiplier_model_z > 0).sum()


# In[15]:


# One off just to get a sense for how many genes are being compared
# Filter genes to only use those shared between our analysis and multiplier
# Check overlap between multiplier genes and our genes
multiplier_genes = list(multiplier_model_z.index)
our_genes = list(pd.read_csv(ls_data_files[0], sep="\t", index_col=0, header=0).index)
shared_genes = set(our_genes).intersection(multiplier_genes)

print(len(our_genes))
print(len(shared_genes))


# In[9]:


# Input: list of generic genes, specific genes, LV matrix
# Compare coverage of generic genes vs LV
# Compare coverage of specific genes vs LV
# Find genes that have a nonzero contribution to LV
# Return number of LV with at least one gene


# In[10]:


# Plot coverage distribution given list of generic coverage, specific coverage
# save plot

