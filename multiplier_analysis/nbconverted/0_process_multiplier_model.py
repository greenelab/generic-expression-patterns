#!/usr/bin/env python
# coding: utf-8

# # Process MULTIPLIER model data
# 
# Raw multiplier model data can be found [here](https://github.com/greenelab/multi-plier).
# 
# This multiplier model (Robject) was formatted in python in [phenoplier repo](https://github.com/greenelab/phenoplier/blob/master/nbs/01_preprocessing/005-multiplier_recount2_models.ipynb). The python pickled data files were output to a [shared google drive](https://drive.google.com/drive/folders/12wvqGzFpFBUOX_CsFkFvvcbKfhl648LE). Since these pickle files were generated using a different conda environment, this notebook is loading the data from the pickle files into .tsv files for use in our analysis.
# 
# **Warning**: This notebook is run using phenoplier environment from [here](https://github.com/greenelab/phenoplier/blob/master/environment/environment.yml)

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import pickle

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


# In[2]:


def read_config(filename):
    """
    Read and parse configuration file containing stored user variables.

    These variables are then passed to the analysis notebooks
    and input to pipeline functions.
    """
    f = open(filename)
    config_dict = {}
    for lines in f:
        items = lines.split("\t", 1)
        config_dict[items[0]] = eval(items[1])
    return config_dict


# In[3]:


# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "human_general_analysis")

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = read_config(config_filename)

local_dir = params["local_dir"]
multiplier_dir = os.path.join(local_dir, "multiplier")


# # Load pickle files and save as .tsv
# 
# The Z is gene x LV matrix. It one of the resulting matrices from The PLIER model, which performs a matrix factorization of gene expression data to get two matrices: loadings (Z: gene x LV) and latent matrix (B: LV x sample). The loadings (Z) are constrained to aligned with curated pathways and gene sets specified by prior knowledge [Figure 1B of Taroni et. al.](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30119-X). This ensure that some but not all latent variables capture known biology. The way PLIER does this is by applying a penalty such that the individual latent variables represent a few gene sets in order to make the latent variables more interpretable. Ideally there would be one latent variable associated with one gene set unambiguously.
# 
# While the PLIER model was trained on specific datasets, MultiPLIER extended this approach to all of recount2, where the latent variables should correspond to specific pathways or gene sets of interest. Therefore, we will look at the coverage of generic genes versus specific genes across these MultiPLIER latent variables in the next [notebook](1_get_multiplier_LV_coverage.ipynb).

# In[4]:


multiplier_z_filename = os.path.join(multiplier_dir, "multiplier_model_z.pkl")


# In[5]:


# Load pickled file
multiplier_model_z = pd.read_pickle(multiplier_z_filename)


# In[7]:


print(multiplier_model_z.shape)
multiplier_model_z.head()


# In[8]:


# make sure I'm seeing the same when loaded with R
assert multiplier_model_z.loc['GAS6', 'LV2'] == 0
assert multiplier_model_z.loc['GAS6', 'LV3'] == 0.039437739697954444
assert multiplier_model_z.loc['GAS6', 'LV984'] == 0.3473620915326928
assert multiplier_model_z.loc['GAS6', 'LV987'] == 0

assert multiplier_model_z.loc['SPARC', 'LV981'] == 0
assert multiplier_model_z.loc['SPARC', 'LV986'].round(8) == 0.12241734


# In[9]:


# Save
multiplier_model_z.to_csv("multiplier_model_z.tsv", sep="\t")


# # Format multiplier summary data
# 
# This summary data matrix contains statistics about each LV - which pathways it was associated with and its significance score. This information is saved in the MultiPLIER model: https://github.com/greenelab/multi-plier/blob/7f4745847b45edf8fef3a49893843d9d40c258cf/23-explore_AAV_recount_LVs.Rmd

# In[10]:


readRDS = ro.r['readRDS']


# In[11]:


multiplier_full_model = readRDS(os.path.join(multiplier_dir,
                                             "recount_PLIER_model.RDS"))


# In[12]:


multiplier_model_matrix = multiplier_full_model.rx2('summary')


# In[13]:


with localconverter(ro.default_converter + pandas2ri.converter):
  multiplier_model_matrix_values = ro.conversion.rpy2py(multiplier_model_matrix)


# In[14]:


multiplier_model_matrix_df = pd.DataFrame(
    data=multiplier_model_matrix_values,
    index=multiplier_model_matrix.rownames,
    columns=multiplier_model_matrix.colnames
)


# In[15]:


multiplier_model_matrix_df.head()


# In[16]:


# Save
multiplier_model_matrix_df.to_csv("multiplier_model_summary.tsv", sep="\t")

