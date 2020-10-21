
# coding: utf-8

# # Process recount2 data
# This notebook does the following:
# 
# 1. Selects template experiment
# 2. Downloads subset of recount2 data, including the template experiment (subset of random experiments + 1 template experiment)
# 3. Train VAE on subset of recount2 data

# In[9]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import numpy as np
import pandas as pd

from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process


# In[10]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_test.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# File that contains gene ranks identified by Crow et. al.
DE_prior_filename = params['reference_gene_filename']

# Template experiment ID
project_id = params['project_id']

# Output file: pickled list of shared genes(generated during gene ID mapping)
shared_genes_filename = params['shared_genes_filename']

# Output files of recount2 template experiment data
raw_template_filename = params['raw_template_filename']
mapped_template_filename = params['mapped_template_filename']
processed_template_filename = params['processed_template_filename']

# Output files of recount2 compendium data
raw_compendium_filename = params['raw_compendium_filename']
mapped_compendium_filename = params['mapped_compendium_filename']
normalized_compendium_filename = params['normalized_compendium_filename']

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params['scaler_filename']


# ## Test: Downloading data

# In[11]:


# Directory where the downloaded files of template experiment will be saved into
template_download_dir = os.path.join(local_dir, "template_download")

# Make sure this directory already exists
os.makedirs(template_download_dir, exist_ok=True)


# In[12]:


get_ipython().run_cell_magic('R', '-i project_id -i template_download_dir -i raw_template_filename -i base_dir', "\nsource(paste(base_dir, 'generic_expression_patterns_modules/download_recount2_data.R', sep='/'))\n\nget_recount2_template_experiment(project_id, template_download_dir, raw_template_filename)")


# ## Test: Renaming gene ids

# In[13]:


# File mapping ensembl ids to hgnc symbols
gene_id_filename = os.path.join(base_dir, dataset_name, "data", "metadata", "ensembl_hgnc_mapping.tsv")


# In[14]:


get_ipython().run_cell_magic('R', '-i raw_template_filename -i gene_id_filename -i base_dir', '\n# Get mapping between ensembl gene ids (ours) to HGNC gene symbols (published)\n# Input: raw_template_filename, output: gene_id_filename\n\nsource(paste(base_dir, \'/generic_expression_patterns_modules/process_names.R\', sep=\'/\'))\n\n# Note: This mapping file from ensembl ids to hgnc symbols is based on the library("biomaRt")\n# that gets updated. In order to get the most up-to-date version, you can delete the \n# ensembl_hgnc_mapping file to re-run the script that generates this mapping.\n\nif (file.exists(gene_id_filename) == FALSE) {\n    get_ensembl_symbol_mapping(raw_template_filename, gene_id_filename)\n}')


# ## Test: processing template data
# This step will map the ensembl gene IDs in raw template data file to hgnc gene symbols, and delete certain columns (genes) and rows (samples). 
# 
# Output files generated in this step: 
# - `shared_genes_filename`: pickled list of shared genes (created only if it doesn't exist yet)
# - `mapped_template_filename`: template data with column names mapped to hgnc gene symbols
# - `processed_template_filename`: template data with some sample rows dropped

# In[7]:


manual_mapping = {                                                                                  
    "ENSG00000187510.7": "PLEKHG7",       
    "ENSG00000230417.11": "LINC00595",                      
    "ENSG00000276085.1": "CCL3L1",                     
    "ENSG00000255374.3": "TAS2R45",                       
}

# metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)

process.process_raw_template_recount2(
    raw_template_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
    mapped_template_filename,
    sample_id_metadata_filename,
    processed_template_filename
)


# In[8]:


# Read data
template_data = pd.read_csv(
    processed_template_filename,
    header=0,
    sep='\t',
    index_col=0
)

# Check samples dropped
print(template_data.shape)
assert(template_data.shape[0] == 24)
template_data.head()


# ## Test: Processing compendium

# In[9]:


process.process_raw_compendium_recount2(
    raw_compendium_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
    mapped_compendium_filename,
    normalized_compendium_filename, 
    scaler_filename
)


# In[10]:


# Check number of genes is equal between the compendium and the template
compendium_data = pd.read_csv(normalized_compendium_filename, sep="\t", index_col=0, header=0)
assert(compendium_data.shape[1] == template_data.shape[1])


# ## Train: VAE training and reproducibility

# In[11]:


# Create VAE directories
output_dirs = [
    os.path.join(base_dir, dataset_name, "models"),
    os.path.join(base_dir, dataset_name, "logs")
]
NN_architecture = params['NN_architecture']

for each_dir in output_dirs:
    new_dir = os.path.join(each_dir, NN_architecture)
    os.makedirs(new_dir, exist_ok=True)


# In[14]:


# Train VAE on new compendium data
train_vae_modules.train_vae(
    config_filename,
    normalized_compendium_filename
)


# In[15]:


# Test reproducibility
expected_log = "data/test_vae_logs.tsv"
actual_log = "logs/NN_2500_30/tybalt_2layer_30latent_stats.tsv"
assert np.all(
    np.isclose(
        pd.read_csv(expected_log, sep="\t").values,
        pd.read_csv(actual_log, sep="\t").values
    )
)

