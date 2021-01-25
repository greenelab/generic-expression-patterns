#!/usr/bin/env python
# coding: utf-8

# # Subset training compendium
# 
# This notebook subsets the normalized compendium to only include those samples that use PAO1 strains. This filtered training compendium will be used to examine the hypothesis that the subset of genes that are found to be generic by SOPHIE are not generic using GAPE-curated experiments because SOPHIE is trained on a compendium containing multiple strains, whereas we suspect that the GAPE experiments are only from a single strain (PAO1). So genes that are generic in other strain contexts will not be detected by GAPE.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd

from ponyo import utils


# ### Set parameters for data processing
# 
# Most parameters are read from `config_filename`. 

# In[2]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_pao1.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
project_id = params['project_id']
raw_compendium_filename = "https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl"

# Load metadata file with annotations per sample
metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "sample_annotations.tsv")

# Load metadata file with annotations per sample
sample_id_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv")

# Output filename
out_compendium_filename = params["raw_compendium_filename"]


# ## Read data and metadata

# In[3]:


# Read normalized compendium
raw_compendium = pd.read_csv(raw_compendium_filename, sep="\t", index_col=0, header=0).T

# Read metadata file
metadata = pd.read_csv(metadata_filename, sep="\t", index_col=0, header=0)

# Read sample id file
project_sample_ids = pd.read_csv(sample_id_filename, sep="\t", index_col=None, header=0)


# In[4]:


raw_compendium.head()


# In[5]:


metadata.head()


# In[6]:


# Get sample ids that have PAO1 strain
pao1_metadata = metadata.loc[metadata["strain"].str.contains("PAO1")]
experiment_ids = list(pao1_metadata.index)

assert "E-GEOD-33245" in experiment_ids


# In[7]:


# Get normalized data associated with PAO1 sample ids
pao1_sample_ids = list(pao1_metadata["ml_data_source"])
print(len(pao1_sample_ids))

pao1_raw_compendium = raw_compendium.loc[pao1_sample_ids]

assert len(pao1_sample_ids) == pao1_raw_compendium.shape[0]


# In[8]:


# Drop samples with no expression activity associated
pao1_raw_compendium.dropna(inplace=True)

print(pao1_raw_compendium.shape)
pao1_raw_compendium


# In[9]:


# Check samples associated with project_id are still in subset
project_sample_ids = list(project_sample_ids["Sample"])

for sample_id in project_sample_ids:
    assert sample_id in list(pao1_raw_compendium.index)


# In[10]:


# Save filtered compendium
pao1_raw_compendium.T.to_csv(out_compendium_filename, sep="\t")

