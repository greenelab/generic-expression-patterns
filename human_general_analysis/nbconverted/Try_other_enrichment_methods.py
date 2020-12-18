#!/usr/bin/env python
# coding: utf-8

# ## Name
# 
# This notebook plugs in other gene set enrichment methods to demonstrate that our method, SOPHIE, can be inserted into different pipelines and work with other methods

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd
import numpy as np
import pickle

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from ponyo import utils
from generic_expression_patterns_modules import calc, process

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
local_dir = params["local_dir"]
project_id = params['project_id']
hallmark_DB_filename = params["pathway_DB_filename"]


# In[4]:


# Load DE stats directory
DE_stats_dir = os.path.join(local_dir, "DE_stats")

# Template experiment DE stats
template_DE_stats_filename = os.path.join(
    DE_stats_dir,
    f"DE_stats_template_data_{project_id}_real.txt"
)


# ## Enrichment methods
# * [ROAST](https://pubmed.ncbi.nlm.nih.gov/20610611/) is available in limma
# * [CAMERA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/) is available in limma
# * [GSVA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/) its own bioconductor package
# * [ORA]() is available in PathwayStudios or David
# 
# TO DO: Write about each method

# In[5]:


# Define function
# ORA works on list of DE
# Apply voom on gene expression >> ROAST, CAMERA, GVSA

# Process data using voom


# Run method on template experiments
# Run method on simulated experiments
# Output table sort by ranking


# In[6]:


# Get summary rank of pathways

