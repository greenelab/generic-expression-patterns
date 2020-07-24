
# coding: utf-8

# # Process pseudomonas data
# This notebook does the following:
# 
# 1. Selects template experiment from the Pseudomonas compendium
# 2. Normalizes the Pseudomonas compendium
# 3. Train VAE on the normalized data

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
import pickle

from ponyo import utils, train_vae_modules, simulate_expression_data
from generic_expression_patterns_modules import process, calc

np.random.seed(123)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

config_file = os.path.abspath(os.path.join(base_dir,
                                           "configs",
                                           "config_pseudomonas_1183.tsv"))
params = utils.read_config(config_file)


# In[3]:


# Load params
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
NN_architecture = params['NN_architecture']
project_id = params['project_id']
metadata_colname = params['metadata_colname']
template_data_file = params['template_data_file']
original_compendium_file = params['compendium_data_file']
normalized_data_file = params['normalized_compendium_data_file']
shared_genes_file = params['shared_genes_file']
scaler_file = params['scaler_transform_file']


# ### Download Pseudomonas compendium
# The compendium is downloaded from https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl

# In[4]:


# Read compendium
original_compendium = pd.read_csv(original_compendium_file,
                                  header=0,
                                  index_col=0,
                                  sep="\t")

if original_compendium.shape != (950, 5549):
    original_compendium = original_compendium.T
    
assert original_compendium.shape == (950, 5549)

print(original_compendium.shape)
original_compendium.head()


# ### Select template experiment
# 
# We manually selected bioproject [E-GEOD-9989](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-9989/?query=George+O%27Toole), which contains 2 samples (3 replicates each) of PA14 WT that are grown on CFBE41o- cells are either treated tobramycin or untreated.
# 
# Another bioproject selected [E-MEXP-1183](https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-1183/), which contains a total of 10 samples. But for now we will select those 4 samples using WT that were measuring the effect of acyl-HSL signal.

# In[5]:


sample_ids = simulate_expression_data.get_sample_ids(project_id, dataset_name, metadata_colname)


# In[6]:


# Get samples from experiment id
template_data = original_compendium.loc[sample_ids]
print(template_data.shape)
template_data.head()


# In[7]:


if project_id == "E-MEXP-1183":
    # drop samples
    sample_ids_to_drop = ["MSC_05.CEL",
                          "MSC_06.CEL",
                          "MSC_07.CEL",
                          "MSC_08.CEL",
                          "MSC_09.CEL",
                          "MSC_10.CEL"]
    template_data = template_data.drop(sample_ids_to_drop)
    
    assert(template_data.shape[0] == 4)


# Note: We are training a compendium using all the samples (including those that are being dropped in the template experiment). However, only the subset of samples (those kept) in the template experiment are those used in the DE analysis in order to ensure the comparison of samples with consistent backgrounds. 
# 
# So there is an inconsistency in the samples used to learn a low-dimensional representation and those used to calculate DE statistics. The inconsistency could possibly effect the DE statistics if the low dimensional space is significantly different including these extra samples vs not. These few samples will likely not effect the space.

# ### Normalize compendium 

# In[8]:


# 0-1 normalize per gene
scaler = preprocessing.MinMaxScaler()
original_data_scaled = scaler.fit_transform(original_compendium)
original_data_scaled_df = pd.DataFrame(original_data_scaled,
                                columns=original_compendium.columns,
                                index=original_compendium.index)

print(original_data_scaled_df.shape)
original_data_scaled_df.head()


# ### Save data files

# In[9]:


# Save data
original_compendium.to_csv(
    original_compendium_file, float_format='%.3f', sep='\t')

template_data.to_csv(template_data_file, float_format='%.5f', sep='\t')

original_data_scaled_df.to_csv(
    normalized_data_file, float_format='%.3f', sep='\t')

# Save scaler transform
outfile = open(scaler_file,'wb')
pickle.dump(scaler,outfile)
outfile.close()

# Save shared genes
# In this case all genes are used
shared_genes = list(original_compendium.columns)

outfile = open(shared_genes_file,'wb')
pickle.dump(shared_genes,outfile)
outfile.close()


# ### Train VAE 

# In[10]:


# Setup directories
# Create VAE directories
output_dirs = [os.path.join(base_dir, dataset_name, "models"),
               os.path.join(base_dir, dataset_name, "logs")]

# Check if analysis output directory exist otherwise create
for each_dir in output_dirs:
    if os.path.exists(each_dir) == False:
        print('creating new directory: {}'.format(each_dir))
        os.makedirs(each_dir, exist_ok=True)

# Check if NN architecture directory exist otherwise create
for each_dir in output_dirs:
    new_dir = os.path.join(each_dir, NN_architecture)
    if os.path.exists(new_dir) == False:
        print('creating new directory: {}'.format(new_dir))
        os.makedirs(new_dir, exist_ok=True)


# In[11]:


# Train VAE on new compendium data
#train_vae_modules.train_vae(config_file,
#                            normalized_data_file)

