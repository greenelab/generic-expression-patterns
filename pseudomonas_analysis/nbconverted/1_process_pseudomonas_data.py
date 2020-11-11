
# coding: utf-8

# # Process pseudomonas data
# This notebook does the following:
# 
# 1. Selects template experiment from the Pseudomonas compendium created from [Tan et. al.](https://msystems.asm.org/content/1/1/e00025-15)
# 2. Normalizes the gene expression data from the Pseudomonas compendium
# 3. Train VAE on the normalized data

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import random
import tensorflow as tf
import pandas as pd
import numpy as np
from sklearn import preprocessing
import pickle

from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process, calc


# In[ ]:


# Set seeds to get reproducible VAE trained models

# The below is necessary in Python 3.2.3 onwards to
# have reproducible behavior for certain hash-based operations.
# See these references for further details:
# https://keras.io/getting-started/faq/#how-can-i-obtain-reproducible-results-using-keras-during-development
# https://docs.python.org/3.4/using/cmdline.html#envvar-PYTHONHASHSEED
# https://github.com/keras-team/keras/issues/2280#issuecomment-306959926

os.environ["PYTHONHASHSEED"] = "0"

# The below is necessary for starting Numpy generated random numbers
# in a well-defined initial state.
np.random.seed(42)

# The below is necessary for starting core Python generated random numbers
# in a well-defined state.
random.seed(12345)

# The below tf.set_random_seed() will make random number generation
# in the TensorFlow backend have a well-defined initial state.
tf.set_random_seed(1234)


# ### Set parameters for data processing
# 
# Most parameters are read from `config_filename`. We manually selected bioproject [GEOD-33245](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-33245/?s_sortby=col_8&s_sortorder=ascending), as the template experiment, which contains multiple different comparisons including WT vs *crc* mutants, WT vs *cbr* mutants in different conditions.

# In[2]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# Column header containing sample ids
metadata_colname = params['metadata_colname']

# Template experiment ID
project_id = params['project_id']

# Output file: pickled list of shared genes(generated during gene ID mapping)
shared_genes_filename = params['shared_genes_filename']

# Output files of pseudomonas template experiment data
raw_template_filename = params['raw_template_filename']
#mapped_template_filename = params['mapped_template_filename']
processed_template_filename = params['processed_template_filename']

# Output files of pseudomonas compendium data
raw_compendium_filename = params['raw_compendium_filename']
processed_compendium_filename = params['processed_compendium_filename']
normalized_compendium_filename = params['normalized_compendium_filename']

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params['scaler_filename']

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv")


# ### Transpose raw pseudomonas compendium and normalize it
# The compendium is from https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl

# In[3]:


process.process_raw_compendium_pseudomonas(
    raw_compendium_filename,
    processed_compendium_filename,
    normalized_compendium_filename,
    scaler_filename,
)


# ### Select template experiment and drop subset of samples
# 
# We manually selected bioproject selected [GEOD-33245](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-33245/?s_sortby=col_8&s_sortorder=ascending), which contains multiple different comparisons including WT vs *crc* mutants, WT vs *cbr* mutants in different conditions.

# In[4]:


process.process_raw_template_pseudomonas(
    processed_compendium_filename,
    project_id,
    dataset_name,
    metadata_colname,
    sample_id_metadata_filename,
    raw_template_filename,
    processed_template_filename,
)


# In[5]:


# Check
template_data = pd.read_csv(
    processed_template_filename, 
    header=0,
    index_col=0,
    sep="\t"
)

if project_id == "E-GEOD-33245":
    assert(template_data.shape[0] == 4)


# **Note:**
# * We are training our VAE model using ALL the data in the compendium.  
# * The template experiment is using a subset of the samples in the real experiment and using those in the DE analysis in order to ensure the comparison of samples with consistent backgrounds (i.e. some experiments have samples with 3 different biological conditions and for now our statistical test is doing a binary comparison).
# * Simulated experiments are generated by shifting the template experiment (using ALL samples in the real experiment) in the latent space. Then dropping the samples to match the template experiment and perform DE analysis.
# 
# 
# So there is an inconsistency in the samples used to learn a low-dimensional representation and those used to calculate DE statistics. This inconsistency should not not change the simulated experiments since all samples in the template experiment are moved the same amount in the latent space. The only way for this inconsistency to effect the simulated experiments is if the low dimensional space is significantly different including all the experiment samples vs only including a subset. However, we believe that such few samples will likely not effect the space. Furthermore, the dataset used to train the VAE should be a general representation of gene expression patterns and shouldn't have to be include the template experiment.

# ### Train VAE 

# In[6]:


# Create VAE directories if needed
output_dirs = [
    os.path.join(base_dir, dataset_name, "models"),
    os.path.join(base_dir, dataset_name, "logs")
]

NN_architecture = params['NN_architecture']

# Check if NN architecture directory exist otherwise create
for each_dir in output_dirs:
    sub_dir = os.path.join(each_dir, NN_architecture)
    os.makedirs(sub_dir, exist_ok=True)


# In[7]:


# Train VAE on new compendium data
train_vae_modules.train_vae(config_filename,
                            normalized_compendium_filename)

