
# coding: utf-8

# # Coverage of MultiPLIER LV
# 
# The goal of this notebook is to examine why genes were found to be generic. Specifically, this notebook is trying to answer the question: Are generic genes found in more multiplier latent variables compared to specific genes?
# 
# The PLIER model performs a matrix factorization of gene expression data to get two matrices: loadings (Z) and latent matrix (B). The loadings (Z) are constrained to aligned with curated pathways and gene sets specified by prior knowledge [Figure 1B of Taroni et. al.](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30119-X). This ensure that some but not all latent variables capture known biology. The way PLIER does this is by applying a penalty such that the individual latent variables represent a few gene sets in order to make the latent variables more interpretable. Ideally there would be one latent variable associated with one gene set unambiguously.
# 
# While the PLIER model was trained on specific datasets, MultiPLIER extended this approach to all of recount2, where the latent variables should correspond to specific pathways or gene sets of interest. Therefore, we will look at the coverage of generic genes versus specific genes across these MultiPLIER latent variables.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import random
import pandas as pd
import seaborn as sns

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ponyo import utils
from generic_expression_patterns_modules import process


# In[2]:


# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "human_general_analysis")

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]


# In[3]:


# Output file
output_figure_filename = "LV_coverage.svg"


# ## Get gene data

# In[4]:


# Get all gene summary file paths
ls_data_files = process.get_gene_summary_files(data_dir)


# In[5]:


# Plot distribution of z-scores to get a sense
data0 = pd.read_csv(ls_data_files[0], sep="\t", index_col=0, header=0)
data1 = pd.read_csv(ls_data_files[1], sep="\t", index_col=0, header=0)
data2 = pd.read_csv(ls_data_files[2], sep="\t", index_col=0, header=0)
data3 = pd.read_csv(ls_data_files[3], sep="\t", index_col=0, header=0)
sns.distplot(data0["abs(Z score)"])
sns.distplot(data1["abs(Z score)"])
sns.distplot(data2["abs(Z score)"])
sns.distplot(data3["abs(Z score)"])


# In[6]:


# Identify generic and specific genes using z-score
# Z-score cutoff was found by calculating invnorm(0.05/17754). 
# To do this in python you can use the following code:
# from scipy.stats import norm
# norm.ppf((0.05/17754)/2)
# Here we are using a p-value = 0.05
# with a Bonferroni correction for 17754 tests, which are
# the number of P. aeruginosa genes

zscore_threshold = 4.68
ls_genes_out = process.get_generic_specific_genes(ls_data_files, zscore_threshold)


# ## Get LV data and filter

# In[7]:


# Load multiplier models
# Converted formatted pickle files (loaded using phenoplier environment) from
# https://github.com/greenelab/phenoplier/blob/master/nbs/01_preprocessing/005-multiplier_recount2_models.ipynb
# into .tsv files
# Raw data was downloaded from https://figshare.com/articles/recount_rpkm_RData/5716033/4
multiplier_model_summary = pd.read_csv("multiplier_model_summary.tsv", sep="\t", index_col=0, header=0)
multiplier_model_z = pd.read_csv("multiplier_model_z.tsv", sep="\t", index_col=0, header=0)


# In[8]:


multiplier_model_summary.head()


# In[9]:


# Only select LVs that are signficantly associated with some pathways or gene set (i.e. FDR < 0.05)
multiplier_model_z_processed = process.process_multiplier_model_z(multiplier_model_z, multiplier_model_summary)


# In[10]:


multiplier_model_z_processed.head()


# ## Quick looks at the data

# In[11]:


# Get a rough sense for how many genes contribute to a given LV
# (i.e. how many genes have a value > 0 for a given LV)
(multiplier_model_z > 0).sum().sort_values(ascending=True)


# In[12]:


# One off just to get a sense for how many genes are being compared
# Filter genes to only use those shared between our analysis and multiplier
# Check overlap between multiplier genes and our genes
multiplier_genes = list(multiplier_model_z.index)
our_genes = list(pd.read_csv(ls_data_files[0], sep="\t", index_col=0, header=0).index)
shared_genes = set(our_genes).intersection(multiplier_genes)

print(len(our_genes))
print(len(shared_genes))


# ## Gene coverage of LV

# In[13]:


generic_cov_ls = []
specific_cov_ls = []

for ifile in range(len(ls_data_files)):
    generic_genes = ls_genes_out[ifile][0]
    specific_genes = ls_genes_out[ifile][1]
    
    # Only include those genes that are in multiplier otherwise will get NAs
    generic_genes_processed, specific_genes_processed = process.process_generic_specific_gene_lists(
        generic_genes, 
        specific_genes, 
        multiplier_model_z_processed
    )
    print(len(generic_genes_processed), len(specific_genes_processed))
    
    generic_cov_i, specific_cov_i = process.get_LV_coverage(
        generic_genes_processed,
        specific_genes_processed,
        multiplier_model_z_processed
    )
    
    generic_cov_ls += list(generic_cov_i)
    specific_cov_ls += list(specific_cov_i)
    
from matplotlib_venn import venn2
venn2([set(generic_cov_ls), set(specific_cov_ls)],
      set_labels=("Generic LVs", "Specific LVs")
     )


# In[14]:


# Create table of unique generic LVs, unique specific LVs, shared LVs
process.create_LV_df(generic_cov_ls, specific_cov_ls, multiplier_model_summary)


# In[15]:


# Find the number of generic and specific genes that have a nonzero contribution to LV
generic_cov = []
specific_cov = []
num_significant_LVs = multiplier_model_z_processed.shape[1]

for ifile in range(len(ls_data_files)):
    generic_genes = ls_genes_out[ifile][0]
    specific_genes = ls_genes_out[ifile][1]
    
    # Only include those genes that are in multiplier otherwise will get NAs
    generic_genes_processed, specific_genes_processed = process.process_generic_specific_gene_lists(
        generic_genes, 
        specific_genes, 
        multiplier_model_z_processed
    )
    print(len(generic_genes_processed), len(specific_genes_processed))
    
    generic_cov_i, specific_cov_i = process.get_LV_coverage(
        generic_genes_processed,
        specific_genes_processed,
        multiplier_model_z_processed
    )
    
    generic_cov.append(len(generic_cov_i)/num_significant_LVs)
    specific_cov.append(len(specific_cov_i)/num_significant_LVs)
    
gene_cov = pd.DataFrame({'Proportion of significantly associated LVs covered': generic_cov + specific_cov,
                         'gene type': ['generic']*len(ls_data_files) + ['specific']*len(ls_data_files) 
                      })


# In[16]:


# Plot coverage distribution given list of generic coverage, specific coverage
print(generic_cov)
print(specific_cov)

import textwrap
fig = sns.swarmplot(data=gene_cov, 
                  x='gene type', 
                  y='Proportion of significantly associated LVs covered', 
                  palette=['grey','powderblue'])
fig.set_xlabel("Gene Type",fontsize=14)
fig.set_ylabel(textwrap.fill("Proportion of pathway-associated LVs covered", width=30),fontsize=14)
fig.tick_params(labelsize=14)
fig.set_title("Coverage of pathway-associated LVs", fontsize=16)


# Its not surprising that there was more coverage for generic genes vs specific genes because there were many more generic genes found. What if we try to randomly sample a generic gene set that is of the same size as the number of specific genes?

# In[17]:


# Randomly generate generic gene list that is the same size as specific gene list
# Re-run coverage analysis
# Find the number of generic and specific genes that have a nonzero contribution to LV
generic_cov_matched = []
specific_cov_matched = []
num_significant_LVs = multiplier_model_z_processed.shape[1]

for ifile in range(len(ls_data_files)):
    generic_genes = ls_genes_out[ifile][0]
    specific_genes = ls_genes_out[ifile][1]
    
    # Only include those genes that are in multiplier otherwise will get NAs
    generic_genes_processed, specific_genes_processed = process.process_generic_specific_gene_lists(
        generic_genes, 
        specific_genes, 
        multiplier_model_z_processed
    )
    
    # Randomly generate generic gene list that is the same size as specific gene list
    generic_genes_processed_subset = random.sample(generic_genes_processed, len(specific_genes_processed))
    print(len(generic_genes_processed_subset), len(specific_genes_processed))
    
    
    generic_cov_i, specific_cov_i = process.get_LV_coverage(
        generic_genes_processed_subset,
        specific_genes_processed,
        multiplier_model_z_processed
    )
    
    generic_cov_matched.append(len(generic_cov_i)/num_significant_LVs)
    specific_cov_matched.append(len(specific_cov_i)/num_significant_LVs)
    
gene_cov_matched = pd.DataFrame({'Proportion of significantly associated LVs covered': generic_cov_matched + specific_cov_matched,
                         'gene type': ['generic']*len(ls_data_files) + ['specific']*len(ls_data_files) 
                      })


# In[18]:


# Plot coverage distribution given list of generic coverage, specific coverage
print(generic_cov_matched)
print(specific_cov_matched)

import textwrap
fig2 = sns.swarmplot(data=gene_cov_matched, 
                  x='gene type', 
                  y='Proportion of significantly associated LVs covered', 
                  palette=['grey','powderblue'])
fig2.set_xlabel("Gene Type",fontsize=14)
fig2.set_ylabel(textwrap.fill("Proportion of pathway-associated LVs covered", width=30),fontsize=14)
fig2.tick_params(labelsize=14)
fig2.set_title("Coverage of pathway-associated LVs", fontsize=16)


# In[19]:


# Save plot
fig.figure.savefig(
        output_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


# **Takeaway:**
# * On average, specific genes cover fewer pathway-associated LVs compared to generic genes, which were found to be linked to all pathway-associated LVs.
# * Warning: This difference in coverage is not surprising given the fact that there 1-6 specific genes identified compared to the ~1000 generic genes found.
# * Some of the LVs that were only found to have generic genes (see [table](generic_only_LV_summary.tsv)) include mainly immune response pathways (monocytes, mast cell activation), wound healing (collagen formation), cell signaling (focal adhesion, integrin1) 
# 
# **Overall, it looks like generic genes are associated with many pathways, acting as *gene hubs*, which is why they are "generic"**
