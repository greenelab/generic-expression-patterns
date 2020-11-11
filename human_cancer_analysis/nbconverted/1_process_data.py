
# coding: utf-8

# # Process Powers et. al. data
# This notebook does the following:
# 
# 1. Select template experiment. This template experiment will be used in the next [notebook](2_identify_generic_genes_pathways.ipynb) to simulate experiments with the same experimental design but testing a different biological process.
# 
# 2.  Powers et. al. data was downloaded from [here](https://www.synapse.org/#!Synapse:syn22685451) and saved to a local directory. Data was processed using steps described below.
# 
# 3. Train VAE on processed data.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import os
import random
import tensorflow as tf
import pandas as pd
import pickle
from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process


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
# Most parameters are read from `config_filename`. We manually selected bioproject [GSE11352](https://www.ncbi.nlm.nih.gov/gds/?term=GSE11352[Accession]) as the template experiment, which contains breast cell lines treated with estradiol at 12H, 24H and 48H.

# In[3]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_cancer.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# File that contains gene ranks identified by Crow et. al.
DE_prior_filename = params['reference_gene_filename']

# Template experiment ID
project_id = params['project_id']

# Column name containing sample ids
metadata_colname = params['metadata_colname']

# Output file: pickled list of shared genes(generated during gene ID mapping)
shared_genes_filename = params['shared_genes_filename']

# Output files of template experiment data
mapped_template_filename = params['mapped_template_filename']
processed_template_filename = params['processed_template_filename']

# Output files of Rani's compendium data
raw_compendium_filename = params['raw_compendium_filename']
mapped_compendium_filename = params['mapped_compendium_filename']
normalized_compendium_filename = params['normalized_compendium_filename']

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params['scaler_filename']


# ### Load compendium data

# In[4]:


raw_compendium = pd.read_csv(raw_compendium_filename, header=0, index_col = 0)
print(raw_compendium.shape)
raw_compendium.head()


# ### Process compendium data
# 
# 1. Drop probe column
# 2. Transpose
# 3. Get only shared genes from Crow et. al.
# 4. Normalize

# In[5]:


# Drop probe column and transpose matrix to be sample x gene
processed_compendium = raw_compendium.drop(columns='Probe').T

# Get only gene expression data for genes in Crow et. al.
our_gene_ids_hgnc = list(processed_compendium.columns)

published_generic_genes = process.get_published_generic_genes(DE_prior_filename)
shared_genes_hgnc = list(
    set(our_gene_ids_hgnc).intersection(published_generic_genes)
)

# In Python, the order of elements in a list that is converted from a set
# is non-deterministic, so it is sorted here to have reproducible result.
shared_genes_hgnc.sort()

# Pickle `shared_genes_hgnc` and save as `shared_genes_filename`
if not os.path.exists(shared_genes_filename):
    with open(shared_genes_filename, 'wb') as pkl_fh:
        pickle.dump(shared_genes_hgnc, pkl_fh)

mapped_compendium = processed_compendium[shared_genes_hgnc]
print(mapped_compendium.shape)
mapped_compendium.head()


# In[6]:


# Save
mapped_compendium.to_csv(mapped_compendium_filename, sep="\t")


# In[7]:


# Normalize data
process.normalize_compendium(
    mapped_compendium_filename,
    normalized_compendium_filename,
    scaler_filename
)


# ### Select and process template data
# 
# 1. Get gene expression associated with `project_id`, which was manually selected by the user and specified in the config file.
# 
# 2. Drop selected samples from template experiments based on metadata, `data/metadata/all_experiments_sample_annotations.csv`, which contains sample comparisons
# 
# Note: The data is not normalized so that we can perform DE analysis in next notebook

# In[8]:


# Note: This is the only notebook using this function, so for now it is included here
# Get sample ids associated with selected project id
def get_sample_ids(experiment_id, mapping_filename):
    """
    Return sample ids for a given experiment id

    """
    # Read in metadata
    metadata = pd.read_csv(mapping_filename, header=0)
    metadata.set_index('gse', inplace=True)
    
    selected_metadata = metadata.loc[experiment_id]
    sample_ids = list(selected_metadata[metadata_colname])

    return sample_ids

# metadata file with mapping from experiment to sample
experiment_to_sample_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "all_experiments_sample_annotations.csv"
)

sample_ids = get_sample_ids(project_id, experiment_to_sample_metadata_filename)

# Get expression data
template_mapped = mapped_compendium.loc[sample_ids]
print(template_mapped.shape)

# Save
template_mapped.to_csv(mapped_template_filename, sep="\t")


# In[9]:


# metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)
    
# Drop sample ids based on metadata file above
sample_ids_to_drop = set()
if os.path.exists(sample_id_metadata_filename):
    # Read in metadata and get samples to be dropped:
    metadata = pd.read_csv(
        sample_id_metadata_filename, sep='\t', header=0, index_col=0
    )
    sample_ids_to_drop = set(metadata[metadata["processing"] == "drop"].index)

# Write the processed recount2 template output file on disk
with open(mapped_template_filename) as ifh, open(processed_template_filename, "w") as ofh:
    for idx, line in enumerate(ifh):
        sample_id = line.split('\t')[0]
        if idx == 0 or sample_id not in sample_ids_to_drop:
            ofh.write(line)


# ### Train VAE 

# In[10]:


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


# In[11]:


# Train VAE on new compendium data
train_vae_modules.train_vae(config_filename, normalized_compendium_filename)

