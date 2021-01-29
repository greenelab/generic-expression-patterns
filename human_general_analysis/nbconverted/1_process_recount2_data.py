#!/usr/bin/env python
# coding: utf-8

# # Process recount2 data
# This notebook does the following:
# 
# 1. Select template experiment. This template experiment will be used in the next [notebook](2_identify_generic_genes_pathways.ipynb) to simulate experiments with the same experimental design but testing a different biological process.
# 
# 
# 2. Download and process SRA data in recount2
#   
#   2a. Download SRA data in recount2 as RangedSummarizedExperiment (rse) object for each project id provided. Raw reads were mapped to genes using Rail-RNA, which includes exon-exon splice junctions. RSE contains counts summarized at the **gene level** using the **Gencode v25 (GRCh38.p7, CHR) annotation** as provided by Gencode.
#   
#   2b. These rse objects return [coverage counts](https://www.bioconductor.org/packages/devel/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html) as   opposed to read counts and therefore we need to apply [scale_counts](https://rdrr.io/bioc/recount/man/scale_counts.html) to scale by **sample coverage** (average number of reads mapped per nucleotide)
#   
#   2c. DESeq performs an internal normalization where geometric mean is calculated for each gene across all samples. The counts for a gene in each sample is then divided by this mean. The median of these ratios in a sample is the size factor for that sample. This procedure corrects for **library size** (i.e. sequencing depth = total number of reads sequenced for a given sample) and RNA composition bias. DESeq expects [un-normalized](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data) data.
# 
# 
# 3. Train VAE on recount2 data

# In[10]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')


# In[11]:


import os
from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process


# In[12]:


# Set seeds to get reproducible VAE trained models
process.set_all_seeds()


# ### Set parameters for data processing
# 
# Most parameters are read from `config_filename`. We manually selected bioproject [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764) as the template experiment, which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients with 2 replicates each.

# In[13]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
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


# ### Download template experiment's expression data and generate raw template data file

# In[14]:


# Directory where the downloaded files of template experiment will be saved into
template_download_dir = os.path.join(local_dir, "template_download")

# Make sure this directory already exists
os.makedirs(template_download_dir, exist_ok=True)


# In[15]:


get_ipython().run_cell_magic('R', '-i project_id -i template_download_dir -i raw_template_filename -i base_dir', "\nsource(paste0(base_dir, '/generic_expression_patterns_modules/download_recount2_data.R'))\n\nget_recount2_template_experiment(project_id, template_download_dir, raw_template_filename)")


# ### Download all recount2 SRA data

# In[7]:


# Directory where the recount2 SRA data files are saved into
compendium_download_dir = os.path.join(local_dir, "compendium_download")
# Make sure this directory already exists
os.makedirs(compendium_download_dir, exist_ok=True)

metadata_dir = local_dir


# In[8]:


get_ipython().run_cell_magic('R', '-i metadata_dir -i compendium_download_dir -i base_dir', "\nsource(paste0(base_dir, '/generic_expression_patterns_modules/download_recount2_data.R'))\n\ndownload_recount2_sra(metadata_dir, compendium_download_dir)")


# ### Create raw recount2 compendium data file
# Compile data in individual projects together into a single raw compendium file. 

# In[9]:


# Output file: `raw_compendium_filename`
process.create_recount2_compendium(compendium_download_dir, raw_compendium_filename)


# ### Subset genes and convert gene names
# For our downstream analysis, we will be comparing our set of differentially expression genes against the set found in [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf), therefore we will limit our genes to include only those genes shared between our starting set of genes and those in publication. 

# In[16]:


# File mapping ensembl ids to hgnc symbols
gene_id_filename = os.path.join(local_dir, "ensembl_hgnc_mapping.tsv")


# In[17]:


get_ipython().run_cell_magic('R', '-i raw_template_filename -i gene_id_filename -i base_dir', '\n# Get mapping between ensembl gene ids (ours) to HGNC gene symbols (published)\n# Input: raw_template_filename, output: gene_id_filename\n\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/process_names.R\'))\n\n# Note: This mapping file from ensembl ids to hgnc symbols is based on the library("biomaRt")\n# that gets updated. In order to get the most up-to-date version, you can delete the \n# ensembl_hgnc_mapping file to re-run the script that generates this mapping.\n\nif (file.exists(gene_id_filename) == FALSE) {\n  get_ensembl_symbol_mapping(raw_template_filename, gene_id_filename)\n}')


# ### Map ensembl gene IDs in template experiment data
# This step will map the ensembl gene IDs in raw template data file to hgnc gene symbols, and delete certain columns (genes) and rows (samples). 
# 
# Output files generated in this step: 
# - `shared_genes_filename`: pickled list of shared genes (created only if it doesn't exist yet)
# - `mapped_template_filename`: template data with column names mapped to hgnc gene symbols

# In[18]:


manual_mapping = {                                                                                  
    "ENSG00000187510.7": "PLEKHG7",       
    "ENSG00000230417.11": "LINC00595",                      
    "ENSG00000276085.1": "CCL3L1",                     
    "ENSG00000255374.3": "TAS2R45",                       
}

process.map_recount2_data(
    raw_template_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    shared_genes_filename,
    mapped_template_filename,
    )


# ### Map ensembl gene IDs in raw compendium file and normalize it
# The mapping process in this step is similar to the one when processing template data. 
# 
# Output files generated in this step:
# - `shared_genes_filename`: pickled list of shared genes (created only if it doesn't exist yet)
# - `mapped_compendium_filename`: compendium data with column names mapped to hgnc gene symbols
# - `normalized_compendium_filename`: normalized compendium data
# - `scaler_filename`: pickled scaler

# In[13]:


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


# ### Train VAE 
# Performed exploratory analysis of compendium data [here](../explore_data/viz_recount2_compendium.ipynb) to help interpret loss curve.

# In[14]:


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


# In[15]:


# Train VAE on new compendium data
train_vae_modules.train_vae(config_filename, normalized_compendium_filename)

