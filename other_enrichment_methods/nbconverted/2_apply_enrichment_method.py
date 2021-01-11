
# coding: utf-8

# ## Apply enrichment method
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
from generic_expression_patterns_modules import ranking

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
statistic = params['gsea_statistic']
hallmark_DB_filename = params["pathway_DB_filename"]
num_runs = params["num_simulated"]
dataset_name = params['dataset_name']

# TO DO:
# What are your choices of methods to use?
enrichment_method = "ROAST"


# In[4]:


# Load DE stats directory
DE_stats_dir = os.path.join(local_dir, "DE_stats")

# Template experiment gene expression
template_expression_filename = os.path.join(base_dir, dataset_name, params["processed_template_filename"])

# Template experiment DE stats
template_DE_stats_filename = os.path.join(
    DE_stats_dir,
    f"DE_stats_template_data_{project_id}_real.txt"
)

# Metadata file with sample grouping to define comparison
metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_groups.tsv"
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
# ORA works on list of DE <<-- how to download and install???

# ROAST, CAMERA <<- what input???

# Process data using voom


# In[6]:


# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "EA_stats"), exist_ok=True)


# In[7]:


# Load pathway data
hallmark_DB_filename = params["pathway_DB_filename"]


# **Apply enrichment to template experiment**
# 
# See supplementary tables: https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz158/5722384

# In[8]:


get_ipython().run_cell_magic('R', '-i base_dir -i local_dir -i project_id -i template_expression_filename -i hallmark_DB_filename -i metadata_filename -i enrichment_method -o template_enriched_pathways', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/other_enrichment_methods.R\'))\n\nout_filename <- paste(local_dir, \n                      "EA_stats/",\n                      enrichment_method,\n                      "_stats_template_data_",\n                      project_id,\n                      "_real.txt", \n                      sep = "")\n\nif (enrichment_method == "GSVA"){\n    \n    template_enriched_pathways <- find_enriched_pathways_GSVA(template_expression_filename, hallmark_DB_filename)\n}\nelse if (enrichment_method == "ROAST"){\n    \n    template_enriched_pathways <- find_enriched_pathways_ROAST(template_expression_filename, metadata_filename, hallmark_DB_filename)\n}\nwrite.table(as.data.frame(template_enriched_pathways), file = out_filename, row.names = F, sep = "\\t")')


# In[9]:


# Quick check
print(template_enriched_pathways.shape)
template_enriched_pathways


# **Apply enrichment to simulated experiments**

# In[10]:


## TO DO: EA stats not outputting in correct location for some reason.


# In[15]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_filename -i metadata_filename -i num_runs -i base_dir -i enrichment_method', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/other_enrichment_methods.R\'))\n\nfor (i in 0:(num_runs-1)){\n    print(i)\n    simulated_expression_filename <- paste(local_dir, \n                                           "pseudo_experiment/selected_simulated_data_",\n                                           project_id,\n                                           "_", \n                                           i,\n                                           "_processed.txt",\n                                           sep = "")\n\n    out_filename <- paste(local_dir,\n                          "EA_stats/",\n                          enrichment_method,\n                          "_stats_simulated_data_",\n                          project_id,\n                          "_",\n                          i,\n                          ".txt", \n                          sep = "")\n    print(out_filename)\n    \n    if (enrichment_method == "GSVA"){\n        enriched_pathways <- find_enriched_pathways_GSVA(simulated_expression_filename, hallmark_DB_filename) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in GSVA")\n    }\n    else if (enrichment_method == "ROAST"){\n        enriched_pathways <- find_enriched_pathways_ROAST(simulated_expression_filename, metadata_filename, hallmark_DB_filename) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in ROAST")\n    }\n}')


# ### TO DO:
# Validate results. Looks like all pathways have the same statistic for ROAST

# ## Format enrichment output
# 
# Each method yields a different output format so we will need to format the data before we can rank and summarize it

# In[12]:


get_ipython().run_cell_magic('R', '-i hallmark_DB_filename -o hallmark_DB_names', 'library("GSA")\n\nhallmark_DB <- GSA.read.gmt(hallmark_DB_filename)\n\nhallmark_DB_names <- as.data.frame(hallmark_DB$geneset.names)')


# In[13]:


ranking.format_enrichment_output(
    local_dir, 
    project_id, 
    enrichment_method, 
    hallmark_DB_names,
    num_runs
)


# ## Rank pathways

# In[14]:


analysis_type = "GSEA"
col_to_rank_pathways = "ES"

template_GSEA_stats_filename = os.path.join(
    local_dir,
    "EA_stats",
    f"{enrichment_method}_stats_template_data_{project_id}_real.txt"    
)
template_GSEA_stats, simulated_GSEA_summary_stats = ranking.process_and_rank_genes_pathways(
    template_GSEA_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_pathways,
)


# ## Pathway summary table

# In[ ]:


# Create intermediate file: "<local_dir>/gene_summary_table_<col_to_rank_pathways>.tsv"
summary_pathway_ranks = ranking.generate_summary_table(
    template_GSEA_stats_filename,
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
    col_to_rank_pathways,
    local_dir,
    'pathway',
    params
)

summary_pathway_ranks.sort_values(by="Z score", ascending=False).head()


# In[ ]:


# Create `pathway_summary_filename`
summary_pathway_ranks.to_csv(pathway_summary_filename, sep='\t')

