
# coding: utf-8

# # Process recount2 data
# This notebook does the following:
# 
# 1. Selects template experiment
# 2. Downloads subset of recount2 data, including the template experiment (subset of random experiments + 1 template experiment)
# 3. Train VAE on subset of recount2 data

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd
import numpy as np
import random
import rpy2
import seaborn as sns
from sklearn import preprocessing
import pickle

from ponyo import generate_template_data, utils, pipeline
from generic_expression_patterns_modules import process, calc

from numpy.random import seed
random_state = 123
seed(random_state)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

config_file = os.path.abspath(os.path.join(base_dir,
                                           "config_human.tsv"))
params = utils.read_config(config_file)


# ### Select template experiment
# 
# We manually selected bioproject [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764), which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients with 2 replicates each.

# In[3]:


# Load params
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
NN_architecture = params['NN_architecture']
project_id = params['project_id']
num_experiments = params['num_experiments']


# ### Download subset of recount2 to use as a compendium
# The compendium will be composed of random experiments + the selected template experiment

# In[4]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Select a\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("derfinder")\n#BiocManager::install("recount")')


# In[5]:


get_ipython().run_cell_magic('R', '', "library('recount')")


# In[6]:


get_ipython().run_cell_magic('R', '-i project_id -i num_experiments -i local_dir -i base_dir', "\nsource('../generic_expression_patterns_modules/download_recount2_data.R')\n\nget_recount2_compendium(project_id, num_experiments, local_dir, base_dir)")


# ### Download expression data for selected project id

# In[7]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir', "\nsource('../generic_expression_patterns_modules/download_recount2_data.R')\n\nget_recount2_template_experiment(project_id, local_dir)")


# ### Subset genes
# For our downstream analysis we will be comparing our set of differentially expression genes against the set found in [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf), therefore we will limit our genes to include only those genes shared between our starting set of genes and those in publication. 

# In[8]:


# Get generic genes identified by Crow et. al.
DE_prior_file = "https://raw.githubusercontent.com/maggiecrow/DEprior/master/DE_Prior.txt"

DE_prior = pd.read_csv(DE_prior_file,
                       header=0,
                       sep="\t")

DE_prior.head()


# In[9]:


# Get list of published generic genes
published_generic_genes = list(DE_prior['Gene_Name'])


# In[10]:


# Get list of our genes
# Load real template experiment
template_data_file = params['template_data_file']

# Read template data
template_data = pd.read_csv(
    template_data_file,
    header=0,
    sep='\t',
    index_col=0)

our_gene_ids = list(template_data.columns)


# In[11]:


# File mapping ensembl ids to hgnc symbols
gene_id_file = os.path.join(
    local_dir,
    "ensembl_hgnc_mapping.tsv")


# In[12]:


get_ipython().run_cell_magic('R', '', 'suppressWarnings(library("biomaRt"))')


# In[13]:


get_ipython().run_cell_magic('R', '-i template_data_file -i gene_id_file', "\n# Get mapping between ensembl gene ids (ours) to HGNC gene symbols (published)\n\nsource('../generic_expression_patterns_modules/process_names.R')\n\nif (file.exists(gene_id_file) == FALSE){\n    gene_id_mapping <- get_ensembl_symbol_mapping(template_data_file, gene_id_file)\n}")


# In[14]:


# Read gene id mapping
gene_id_mapping = pd.read_csv(
        gene_id_file,
        header=0,
        sep='\t',
        index_col=0)

print(gene_id_mapping.shape)
gene_id_mapping.head()


# In[15]:


# Get mapping between ensembl ids with and without version numbers
# Expressiond data uses ensembl gene ids with version number 
ensembl_gene_ids = pd.DataFrame(data={'ensembl_version': our_gene_ids,
                                      'ensembl_parsed': [gene_id.split('.')[0] for gene_id in our_gene_ids]})

print(ensembl_gene_ids.shape)
ensembl_gene_ids.head()


# In[16]:


# Map ensembl ids with version number to gene_id_mapping_df
gene_id_mapping = pd.merge(gene_id_mapping, 
                           ensembl_gene_ids, 
                           left_on='ensembl_gene_id',
                           right_on='ensembl_parsed', 
                           how='outer')

print(gene_id_mapping.shape)
gene_id_mapping.set_index('ensembl_version', inplace=True)
gene_id_mapping.head()


# Since this experiment contains both RNA-seq and smRNA-seq samples which are in different ranges so we will drop smRNA samples so that samples are within the same range. The analysis identifying these two subsets of samples can be found in this [notebook](0_explore_input_data.ipynb)

# In[17]:


# Replace ensembl ids with gene symbols
template_data = process.replace_ensembl_ids(template_data,
                                            gene_id_mapping)


# In[18]:


template_data.head()


# In[19]:


# Get intersection of gene lists
our_gene_ids_hgnc = template_data.columns
shared_genes_hgnc = list(set(our_gene_ids_hgnc).intersection(published_generic_genes))
print(len(shared_genes_hgnc))


# In[20]:


# Save shared genes
shared_genes_file = os.path.join(
    local_dir,
    "shared_gene_ids.pickle")

outfile = open(shared_genes_file,'wb')
pickle.dump(shared_genes_hgnc,outfile)
outfile.close()


# In[21]:


# Drop smRNA samples so that samples are within the same range
smRNA_samples = ["SRR493961",
                 "SRR493962",
                 "SRR493963",
                 "SRR493964",
                 "SRR493965",
                 "SRR493966",
                 "SRR493967",
                 "SRR493968",
                 "SRR493969",
                 "SRR493970",
                 "SRR493971",
                 "SRR493972"]


# In[22]:


# Drop samples
template_data = template_data.drop(smRNA_samples)


# In[23]:


# Drop genes
template_data = template_data[shared_genes_hgnc]

print(template_data.shape)
template_data.head()


# In[24]:


print(len(template_data.columns) - len(shared_genes_hgnc))


# *Note:* There is a difference in the number of `shared_genes_hgnc` and genes in the template experiment because 3 genes have 2 different ensembl gene ids have map to the same hgnc symbol (one forward, one reverse)

# In[25]:


# Save 
template_data.to_csv(template_data_file, float_format='%.5f', sep='\t')


# ### Normalize compendium 

# In[26]:


# Load real gene expression data
original_compendium_file = params['compendium_data_file']


# In[27]:


# Read data
original_compendium = pd.read_table(
    original_compendium_file,
    header=0,
    sep='\t',
    index_col=0)

print(original_compendium.shape)
original_compendium.head()


# In[28]:


# Replace ensembl ids with gene symbols
original_compendium = process.replace_ensembl_ids(original_compendium,
                                                gene_id_mapping)


# In[29]:


# Drop genes
original_compendium = original_compendium[shared_genes_hgnc]

original_compendium.head()


# In[30]:


# 0-1 normalize per gene
scaler = preprocessing.MinMaxScaler()
original_data_scaled = scaler.fit_transform(original_compendium)
original_data_scaled_df = pd.DataFrame(original_data_scaled,
                                columns=original_compendium.columns,
                                index=original_compendium.index)

print(original_data_scaled_df.shape)
original_data_scaled_df.head()


# In[31]:


# Save data
normalized_data_file = params['normalized_compendium_data_file']

original_data_scaled_df.to_csv(
    normalized_data_file, float_format='%.3f', sep='\t')

original_compendium.to_csv(
    original_compendium_file, float_format='%.3f', sep='\t')

# Save scaler transform
scaler_file = params['scaler_transform_file']

outfile = open(scaler_file,'wb')
pickle.dump(scaler,outfile)
outfile.close()


# ### Train VAE 
# Performed exploratory analysis of compendium data [here](../explore_data/viz_recount2_compendium.ipynb) to help interpret loss curve.

# In[32]:


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


# In[33]:


# Train VAE on new compendium data
pipeline.train_vae(config_file,
                   normalized_data_file)

