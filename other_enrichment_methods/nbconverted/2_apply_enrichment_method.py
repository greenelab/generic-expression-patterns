
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
project_id = params["project_id"]
statistic = params["gsea_statistic"]
hallmark_DB_filename = params["pathway_DB_filename"]
num_runs = params["num_simulated"]
dataset_name = params["dataset_name"]

# Select enrichment method
# enrichment_method = ["GSEA", GSVA", "ROAST", "CAMERA", "ORA"]
# If enrichment_method == "GSEA" then use "padj" to rank
# If enrichment_method == "GSVA" then use "ES" to rank
# If enrichment_method == "ROAST" or "CAMERA" then use "FDR" to rank
# If using "ORA" then use "p.adjust" to rank
enrichment_method = "ROAST"
col_to_rank_pathways = params["rank_pathways_by"]


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


# In[5]:


# Output file
pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}_{enrichment_method}.tsv"
)


# ## Enrichment methods
# These methods are important because they permit differential expression questions to be posed in terms of ensembles of genes representing pathways or other biologically interpretable processes. 
# 
# There are two main groups of methods: ‘self-contained’ gene set tests examine a set of genes in their own right without reference to other genes in the genome whereas ‘competitive’ gene set tests compare genes in the test set relative to all other genes. Self-contained tests are of interest for assessing the relevance of an individual biological process to the experiment at hand, whereas the competitive tests focus more on distinguishing the most important biological processes from those that are less important. Competitive tests are overwhelmingly more commonly used in the genomic literature. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/)
# 
# * [ROAST](https://pubmed.ncbi.nlm.nih.gov/20610611/) (rotation gene set tests) performs a focused gene set testing, in which interest focuses on a few gene sets as opposed to a large dataset. (available in limma).
#   * Self contained gene set test
#   * Instead of permutations they use rotations (i.e. fractional permutation) in order to allow for more complex experimental designs than binary experiments (i.e. time-course, more than 2 groups)
#   * Also esimates correlation between genes
#   * Best power to detect modest fold changes with minority of genes changed
# * [CAMERA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/) (Correlation Adjusted MEan RAnk gene set test) is based on the idea of estimating the variance inflation factor associated with inter-gene correlation, and incorporating this into parametric or rank-based test procedures. (available in limma) 
#   * Competitive gene set test
#   * Performs the same rank-based test procedure as GSEA, but also estimates the correlation between genes, instead of treating genes as independent
#   * Recall GSEA: 1) Rank all genes using DE association statistics. 2) An enrichment score (ES) is defined as the maximum distance from the middle of the ranked list. Thus, the enrichment score indicates whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list (indicating a correlation with change in expression). 3) Estimate the statistical significance of the ES by a phenotypic-based permutation (permute samples assigned to label) test in order to produce a null distribution for the ES (i.e. scores based on permuted phenotype)
#   * Appropriate for small and large fold changes
# * [GSVA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/) (Gene Set Variation Analysis) calculates sample-wise gene set enrichment scores as a function of genes inside and outside the gene set. This method is well-suited for assessing gene set variation across a dichotomous phenotype. (biocontuctor package GSVA) 
#   * Competitive gene set test
#   * Estimates variation of gene set enrichment over the samples independently of any class label
# * [ORA](https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enricher) (over-representation analysis) uses the hypergeometric test to determine if there a significant over-representation of pathway in the selected set of DEGs. Here we're using clusterProfiler library but there are multiple options for this analysis. See [slide 6](https://docs.google.com/presentation/d/1t4rK7UiLAeIKIzeRJK-YzspNUfGM-8nuRCcevh2lx34/edit?usp=sharing)

# In[6]:


# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "GSA_stats"), exist_ok=True)


# In[7]:


# Load pathway data
hallmark_DB_filename = params["pathway_DB_filename"]


# **Apply enrichment to template experiment**
# 
# See supplementary tables: https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz158/5722384

# In[8]:


get_ipython().run_cell_magic('R', '-i base_dir -i local_dir -i project_id -i template_expression_filename -i hallmark_DB_filename -i metadata_filename -i enrichment_method -o template_enriched_pathways', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/other_enrichment_methods.R\'))\n\nout_filename <- paste(local_dir, \n                      "GSA_stats/",\n                      enrichment_method,\n                      "_stats_template_data_",\n                      project_id,\n                      "_real.txt", \n                      sep = "")\n\nif (enrichment_method == "GSVA"){\n    \n    template_enriched_pathways <- find_enriched_pathways_GSVA(\n        template_expression_filename,\n        hallmark_DB_filename\n    )\n}\nelse if (enrichment_method == "ROAST"){\n    \n    template_enriched_pathways <- find_enriched_pathways_ROAST(\n        template_expression_filename,\n        metadata_filename,\n        hallmark_DB_filename\n    )\n}\nelse if (enrichment_method == "CAMERA"){\n    \n    template_enriched_pathways <- find_enriched_pathways_CAMERA(\n        template_expression_filename,\n        metadata_filename,\n        hallmark_DB_filename\n    )\n}\nelse if (enrichment_method == "ORA"){\n    \n    template_enriched_pathways <- find_enriched_pathways_ORA(\n        template_expression_filename,\n        metadata_filename, \n        hallmark_DB_filename\n    )\n}\nwrite.table(as.data.frame(template_enriched_pathways), file = out_filename, row.names = F, sep = "\\t")')


# In[9]:


# Quick check
print(template_enriched_pathways.shape)
template_enriched_pathways


# **Apply enrichment to simulated experiments**
# 
# Note: GSA takes a while to run (each simulated experiment took ~2-3 minutes)

# In[10]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_filename -i metadata_filename -i num_runs -i base_dir -i enrichment_method', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/other_enrichment_methods.R\'))\n\nfor (i in 0:(num_runs-1)){\n    simulated_expression_filename <- paste(local_dir, \n                                           "pseudo_experiment/selected_simulated_data_",\n                                           project_id,\n                                           "_", \n                                           i,\n                                           "_processed.txt",\n                                           sep = "")\n\n    out_filename <- paste(local_dir,\n                          "GSA_stats/",\n                          enrichment_method,\n                          "_stats_simulated_data_",\n                          project_id,\n                          "_",\n                          i,\n                          ".txt", \n                          sep = "")\n    \n    if (enrichment_method == "GSVA"){\n        enriched_pathways <- find_enriched_pathways_GSVA(\n            simulated_expression_filename, \n            hallmark_DB_filename\n        ) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in GSVA")\n    }\n    else if (enrichment_method == "ROAST"){\n        enriched_pathways <- find_enriched_pathways_ROAST(\n            simulated_expression_filename,\n            metadata_filename,\n            hallmark_DB_filename\n        ) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in ROAST")\n    }\n    else if (enrichment_method == "CAMERA"){\n        enriched_pathways <- find_enriched_pathways_CAMERA(\n            simulated_expression_filename,\n            metadata_filename, \n            hallmark_DB_filename\n        ) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in CAMERA")\n    }\n    else if (enrichment_method == "ORA"){\n        enriched_pathways <- find_enriched_pathways_ORA(\n            simulated_expression_filename,\n            metadata_filename, \n            hallmark_DB_filename\n        ) \n        write.table(as.data.frame(enriched_pathways), file = out_filename, row.names = F, sep = "\\t")\n        print("in ORA")\n    }\n}')


# ## Format enrichment output
# 
# Each method yields a different output format so we will need to format the data before we can rank and summarize it

# In[11]:


get_ipython().run_cell_magic('R', '-i hallmark_DB_filename -o hallmark_DB_names', 'library("GSA")\n\nhallmark_DB <- GSA.read.gmt(hallmark_DB_filename)\n\nhallmark_DB_names <- as.data.frame(hallmark_DB$geneset.names)')


# In[12]:


ranking.format_enrichment_output(
    local_dir, 
    project_id, 
    enrichment_method, 
    hallmark_DB_names,
    num_runs
)


# ## Rank pathways

# In[13]:


analysis_type = "GSA"

template_GSEA_stats_filename = os.path.join(
    local_dir,
    "GSA_stats",
    f"{enrichment_method}_stats_template_data_{project_id}_real.txt"    
)
template_GSEA_stats, simulated_GSEA_summary_stats = ranking.process_and_rank_genes_pathways(
    template_GSEA_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_pathways,
    enrichment_method
)


# ## Pathway summary table

# In[14]:


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


# In[15]:


# Create `pathway_summary_filename`
summary_pathway_ranks.to_csv(pathway_summary_filename, sep='\t')

