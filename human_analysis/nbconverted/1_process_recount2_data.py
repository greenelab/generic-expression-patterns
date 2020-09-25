
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

# In[ ]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')


# In[ ]:


import os
import sys
import pandas as pd
import numpy as np
import pickle

from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process


# ### Set parameters for data processing
# 
# Most parameters are read from `config_file`. We manually selected bioproject [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764) as the template experiment, which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients with 2 replicates each.

# In[ ]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_file = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human.tsv")
)

params = utils.read_config(config_file)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# File that contains generic genes identified by Crow et. al.
DE_prior_filename = params['reference_gene_file']

# Template experiment ID
project_id = params['project_id']

# Output files for recount2 template experiment
raw_template_filename = params['raw_template_filename']
mapped_template_filename = params['mapped_template_filename']
processed_template_filename = params['processed_template_filename']

# Output files for recount2 compendium data
raw_compendium_filename = params['raw_compendium_filename']
mapped_compendium_filename = params['mapped_compendium_filename']
normalized_compendium_filename = params['normalized_compendium_filename']

# Are these two pickle files for debugging only?
shared_genes_file = params['shared_genes_file']
scaler_file = params['scaler_transform_file']


# ### Download template experiment's expression data and generate raw template data file

# In[ ]:


# Directory where the downloaded files of template experiment will be saved into
template_download_dir = os.path.join(local_dir, "template_download")

# Make sure this directory already exists
os.makedirs(template_download_dir, exist_ok=True)


# In[ ]:


get_ipython().run_cell_magic('R', '-i project_id -i template_download_dir -i raw_template_filename', "\nsource('../generic_expression_patterns_modules/download_recount2_data.R')\n\nget_recount2_template_experiment(project_id, template_download_dir, raw_template_filename)")


# ### Download all recount2 SRA data

# In[ ]:


# Directory where the recount2 SRA data files are saved into
compendium_download_dir = os.path.join(local_dir, "compendium_download")
# Make sure this directory already exists
os.makedirs(compendium_download_dir, exist_ok=True)

metadata_dir = local_dir


# In[ ]:


get_ipython().run_cell_magic('R', '-i metadata_dir -i compendium_download_dir', "\nsource('../generic_expression_patterns_modules/download_recount2_data.R')\n\ndownload_recount2_sra(metadata_dir, compendium_download_dir)")


# ### Create raw recount2 compendium data file
# Compile data in individual projects together into a single raw compendium file. 

# In[ ]:


# Output file: `raw_compendium_filename`
process.create_recount2_compendium(compendium_download_dir, raw_compendium_filename)


# ### Subset genes and convert gene names
# For our downstream analysis, we will be comparing our set of differentially expression genes against the set found in [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf), therefore we will limit our genes to include only those genes shared between our starting set of genes and those in publication. 

# In[ ]:


# File mapping ensembl ids to hgnc symbols
gene_id_filename = os.path.join(local_dir, "ensembl_hgnc_mapping.tsv")


# In[ ]:


get_ipython().run_cell_magic('R', '', 'suppressWarnings(library("biomaRt"))')


# In[ ]:


get_ipython().run_cell_magic('R', '-i raw_template_filename -i gene_id_filename', '\n# Get mapping between ensembl gene ids (ours) to HGNC gene symbols (published)\n# Input: raw_template_filename, output: gene_id_filename\n\nsource(\'../generic_expression_patterns_modules/process_names.R\')\n\n# Note: This mapping file from ensembl ids to hgnc symbols is based on the library("biomaRt")\n# that gets updated. In order to get the most up-to-date version, you can delete the \n# ensembl_hgnc_mapping file to re-run the script that generates this mapping.\n\nif (file.exists(gene_id_filename) == FALSE) {\n  get_ensembl_symbol_mapping(raw_template_filename, gene_id_filename)\n}')


# ### Map ensembl gene IDs in template experiment data
# This step will map the ensembl gene IDs in raw template data file to hgnc gene symbols, and delete certain columns (genes) and rows (samples). 
# 
# Two output files will be generated in this step: `mapped_template_filename` and `processed_template_filename`.

# In[ ]:


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

process.process_raw_template(
    raw_template_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    mapped_template_filename,
    sample_id_metadata_filename,
    processed_template_filename
)


# ### Map ensembl gene IDs in raw compendium file and normalize it
# The mapping process in this step is similar to the one when processing template data. 
# 
# Two output files will be generated in this step: `mapped_compendium_filename` and `normalized_compendium_filename`.

# In[ ]:


process.process_raw_compendium(
    raw_compendium_filename,
    gene_id_filename,
    manual_mapping,
    DE_prior_filename,
    mapped_compendium_filename,
    normalized_compendium_filename
)


# ### Train VAE 
# Performed exploratory analysis of compendium data [here](../explore_data/viz_recount2_compendium.ipynb) to help interpret loss curve.

# In[ ]:


# Create VAE directories if needed
output_dirs = [
    os.path.join(base_dir, dataset_name, "models"),
    os.path.join(base_dir, dataset_name, "logs")
]

NN_architecture = params['NN_architecture']

# Check if NN architecture directory exist otherwise create
for each_dir in output_dirs:
    new_dir = os.path.join(each_dir, NN_architecture)
    os.makedirs(new_dir, exist_ok=True)


# In[ ]:


# Train VAE on new compendium data
train_vae_modules.train_vae(config_file, normalized_compendium_filename)

