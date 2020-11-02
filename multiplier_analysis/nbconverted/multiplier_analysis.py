
# coding: utf-8

# # Multiplier analysis
# 
# The goal of this notebook is to examine why genes were found to be generic. Specifically, this notebook is trying to answer the question: Are generic genes found in more multiplier latent variables compared to specific genes?
# 
# The PLIER model performs a matrix factorization of gene expression data to get two matrices: loadings (Z) and latent matrix (B). The loadings (Z) are constrained to aligned with curated pathways and gene sets specified by prior knowledge [Figure 1B of Traoni et. al.](). This ensure that some but not all latent variables capture known biology. The way PLIER does this is by applying a penalty such that the individual latent variables represent a few gene sets in order to make the latent variables more interpretable. Ideally there would be one latent variable associated with one gene set unambiguously.
# 
# While the PLIER model was trained on specific datasets, MULTIPLIER extended this approach to all of recount2, where the latent variables should correspond to specific pathways or gene sets of interest. Therefore, we will look at the coverage of generic genes versus specific genes across these MULTIPLIER latent variables.

# In[19]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import seaborn as sns

from generic_expression_patterns_modules import process


# In[2]:


# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "human_general_analysis")


# In[ ]:


# Output file
output_figure_filename = "LV_coverage.svg"


# ## Get gene data

# In[3]:


# Get all gene summary file paths
ls_data_files = process.get_gene_summary_files(data_dir)


# In[5]:


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
ls_genes_out = process.get_generic_specific_genes(ls_data_files, zscore_threshold)

# TO DO:
# Is this how we want to define the genes? ranking?


# ## Get LV data and filter

# In[ ]:


# TO DO
# Filter multiplier matrix somehow


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


# ## Quick looks at the data

# In[8]:


# Get a rough sense for how many genes contribute to a given LV
# (i.e. how many genes have a value > 0 for a given LV)
(multiplier_model_z > 0).sum()


# In[9]:


# Get a rough sense for how many LV have a nonzero association with a known pathway/gene set
# (i.e. how many genes have a value > 0 per LV)
((multiplier_model_u > 0).sum()>0).sum()


# In[10]:


# One off just to get a sense for how many genes are being compared
# Filter genes to only use those shared between our analysis and multiplier
# Check overlap between multiplier genes and our genes
multiplier_genes = list(multiplier_model_z.index)
our_genes = list(pd.read_csv(ls_data_files[0], sep="\t", index_col=0, header=0).index)
shared_genes = set(our_genes).intersection(multiplier_genes)

print(len(our_genes))
print(len(shared_genes))


# ## Gene coverage of LV

# In[17]:


# Find the number of generic and specific genes that have a nonzero contribution to LV
generic_cov = []
specific_cov = []

for ifile in range(len(ls_data_files)):
    generic_genes = ls_genes_out[ifile][0]
    specific_genes = ls_genes_out[ifile][1]
    
    generic_cov_i, specific_cov_i = process.get_LV_coverage(generic_genes, specific_genes, multiplier_model_z)
    
    generic_cov.append(generic_cov_i)
    specific_cov.append(specific_cov_i)
    
gene_cov = pd.DataFrame({'coverage': generic_cov + specific_cov,
                         'gene type': ['generic']*len(ls_data_files) + ['specific']*len(ls_data_files) 
                      })


# In[21]:


# Plot coverage distribution given list of generic coverage, specific coverage
print(generic_cov)
print(specific_cov)

sns.boxplot(data=gene_cov, x='gene type', y='coverage')


# In[ ]:


# Save plot
fig.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )

