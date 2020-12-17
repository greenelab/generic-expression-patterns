#!/usr/bin/env python
# coding: utf-8

# # Debug
# 
# During a call with Casey and Jim, they noticed 2 unusual things in the generic_gene_summary table:
# * Not all values in the `num_simulated` column were equal to 25, which should be the case
# * There are some genes that do not have any DE statistics reported. This is the case using a template experiment that is included in the recount2 training compendium, so its not just an issue using an external template experiment.
# 
# These two issues were NOT observed in the other generic_gene_summary tables training compendium = Powers et. al. or Pseudomonas datasets.
# 
# See [example summary tables](https://docs.google.com/spreadsheets/d/1aqSPTLd5bXjYOBoAjG7bM8jpy5ynnCpG0oBC2zGzozc/edit?usp=sharing)
# 
# Given that this is only observed in the recount2 training dataset, we suspect that this is an issue with mapping genes from ensembl ids (raw data) to hgnc ids (needed to compare against validation dataset)

# In[35]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import re
import pandas as pd
from ponyo import utils
from generic_expression_patterns_modules import process, new_experiment_process


# In[2]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)


# In[3]:


# Load params
# Output files of recount2 template experiment data
raw_template_recount2_filename = params['raw_template_filename']
mapped_template_recount2_filename = os.path.join(
    base_dir, 
    "human_general_analysis",
    params['processed_template_filename']
)

# Local directory to store intermediate files
local_dir = params['local_dir']

# ID for template experiment
# This ID will be used to label new simulated experiments
project_id = params['project_id']


# In[4]:


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_new_experiment.tsv")
)

params = utils.read_config(config_filename)


# In[5]:


# Load params
# Output files of recount2 template experiment data
raw_template_nonrecount2_filename = params['raw_template_filename']
mapped_template_nonrecount2_filename = params['processed_template_filename']


# ## Read data

# In[6]:


# Read raw template data
# This is the data **before** gene ids were mapped
raw_template_recount2 = pd.read_csv(raw_template_recount2_filename, sep="\t", index_col=0, header=0)
raw_template_nonrecount2 = pd.read_csv(raw_template_nonrecount2_filename, sep="\t", index_col=0, header=0).T

# Read mapped template data
# This is the data **after** gene ids were mapped
mapped_template_recount2 = pd.read_csv(mapped_template_recount2_filename, sep="\t", index_col=0, header=0)
mapped_template_nonrecount2 = pd.read_csv(mapped_template_nonrecount2_filename, sep="\t", index_col=0, header=0)


# In[7]:


print(raw_template_recount2.shape)
raw_template_recount2.head()


# In[8]:


print(mapped_template_recount2.shape)
mapped_template_recount2.head()


# In[9]:


print(raw_template_nonrecount2.shape)
raw_template_nonrecount2.head()


# In[10]:


print(mapped_template_nonrecount2.shape)
mapped_template_nonrecount2.head()


# ## Look up some genes

# ## Case 1:
# Genes have only simulated statistics using recount2 template experiment:
# * ENSG00000169717 --> ACTRT2
# * ENSG00000184895 --> SRY

# In[11]:


ensembl_ids = ["ENSG00000169717",
               "ENSG00000184895",
               "ENSG00000124232",
               "ENSG00000261713",
               "ENSG00000186818",
               "ENSG00000160882"               
              ]

ensembl_version_ids = raw_template_recount2.columns

for sub in ensembl_ids:
    for x in ensembl_version_ids:
        if re.search(sub, x):
            print(x)


# In[12]:


raw_template_recount2[["ENSG00000169717.6", "ENSG00000184895.7"]].sum()


# In[13]:


mapped_template_recount2[["ACTRT2","SRY"]].sum()


# Looks like reason for missing values in template experiment could be due to having all 0 counts

# ## Case 2:
# 
# Genes that have all statistics present but number of simulated experiments < 25 using recount2 template. These genes are also only have simulated statistics using non-recount2 template experiment:
# * ENSG00000186818 --> LILRB4
# * ENSG00000160882 --> CYP11B1

# In[14]:


raw_template_nonrecount2[["LILRB4","CYP11B1"]].sum()


# In[15]:


mapped_template_nonrecount2[["LILRB4","CYP11B1"]].sum()


# So far, it seems that those genes that are missing template statistics have all 0 counts in the template experiment.

# In[16]:


raw_template_recount2[["ENSG00000186818.12",
                       "ENSG00000160882.11"]].sum()


# In[17]:


mapped_template_recount2[["LILRB4","CYP11B1"]].sum()


# Overall there isn't a trend found in these genes missing some number of simulated experiments, so let's try looking at the simulated experiments. At this point we suspect that the missing simulated experiments are those where genes have all 0 counts.

# In[18]:


# Get list of files
simulated_dir = os.path.join(
    local_dir,
    "pseudo_experiment",
)

simulated_filename_list = []
for file in os.listdir(simulated_dir):
    if (project_id in file) and ("simulated" in file) and ("encoded" not in file):
        simulated_filename_list.append(os.path.join(simulated_dir,file))


# In[19]:


assert len(simulated_filename_list) ==25


# In[20]:


# For each simulated experiment, check how many have all 0 counts for this gene
# Is this number the same number of missing simulated experiments?
counter_LILRB4 = 0
counter_CYP11B1 = 0
for filename in simulated_filename_list:
    simulated_data = pd.read_csv(filename, sep="\t", index_col=0, header=0)
    if simulated_data["LILRB4"].sum() == 0:
        counter_LILRB4 += 1
    if simulated_data["CYP11B1"].sum() == 0:
        counter_CYP11B1 += 1
        
# Verified LILRB4 to be missing 2 experiments (23 total experiments)
# Verified CYP11B1 to be missing 8 experiments (17 total experiments)
# Can look this up in google sheet
print(counter_LILRB4, counter_CYP11B1)


# ## Case 3:
# 
# Genes that do not have a p-value using non-recount2 template experiment:
# * ENSG00000124232 --> RBPJL
# * ENSG00000261713 --> SSTR5-AS1

# Following the theme observed in Case 1 and 2, we suspect that these "missing" p-values actually indicate p-value =0

# In[21]:


# Look up these values in DE statis output files
DE_stats_dir = os.path.join(
    local_dir,
    "DE_stats"
)

template_nonrecount2_DE_stats_filename = os.path.join(
    DE_stats_dir,
    "DE_stats_template_data_cis-par-KU1919_real.txt"
)

template_nonrecount2_DE_stats = pd.read_csv(
    template_nonrecount2_DE_stats_filename,
    sep="\t",
    index_col=0,
    header=0
)

template_nonrecount2_DE_stats.loc["RBPJL"]


# In[22]:


template_nonrecount2_DE_stats.loc["SSTR5.AS1"]


# According to this [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA), genes with NA adjusted p-values indicate those genes that have been automatically filtered by DESeq2 as very likely to not have no significance. 
# 
# Let's test calling DESeq2, turning off filtering

# In[23]:


transposed_template_filename = "/home/alexandra/Documents/Data/Generic_expression_patterns/Costello_BladderCancer_ResistantCells_Counts_12-8-20_transposed.txt"

new_experiment_process.transpose_save(raw_template_nonrecount2_filename, transposed_template_filename)


# In[29]:


project_id = "cis-par-KU1919"
local_dir = "/home/alexandra/Documents/Data/Generic_expression_patterns/"
mapped_compendium_filename = "/home/alexandra/Documents/Data/Generic_expression_patterns/mapped_recount2_compendium.tsv"


# In[27]:


# Check that the feature space matches between template experiment and VAE model.  
# (i.e. ensure genes in template and VAE model are the same).
mapped_template_experiment = new_experiment_process.compare_match_features(
    transposed_template_filename,
    mapped_compendium_filename
)
mapped_template_filename = transposed_template_filename


# In[30]:


# Load metadata file with processing information
sample_id_metadata_filename = os.path.join(
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)

# Read in metadata
metadata = pd.read_csv(sample_id_metadata_filename, sep='\t', header=0, index_col=0)

# Get samples to be dropped
sample_ids_to_drop = list(metadata[metadata["processing"] == "drop"].index)


# In[31]:


# Modify template experiment
process.subset_samples_template(
    mapped_template_filename,
    sample_ids_to_drop,
)


# In[32]:


process.recast_int_template(mapped_template_filename)


# In[33]:


# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    "data",
    "metadata",
    f"{project_id}_groups.tsv"
)

# Check whether ordering of sample ids is consistent between gene expression data and metadata
process.compare_and_reorder_samples(mapped_template_filename, metadata_filename)


# In[39]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i mapped_template_filename -i local_dir -i base_dir', 'library("limma")\nlibrary("DESeq2")\n# Manually change DESeq2 call with additional parameter\n\nget_DE_stats_DESeq <- function(metadata_file,\n                               experiment_id,\n                               expression_file,\n                               data_type,\n                               local_dir,\n                               run) {\n\n  # This function performs DE analysis using DESeq.\n  # Expression data in expression_file are grouped based on metadata_file\n  #\n  # Arguments\n  # ---------\n  # metadata_file: str\n  #   File containing mapping between sample id and group\n  #\n  # experiment_id: str\n  #   Experiment id used to label saved output filee\n  #\n  # expression_file: str\n  #   File containing gene expression data\n  #\n  # data_type: str\n  #   Either \'template\' or \'simulated\' to label saved output file\n  #\n  # local_dir: str\n  #   Directory to save output files to\n  #\n  # run: str\n  #   Used as identifier for different simulated experiments\n\n  expression_data <- t(as.matrix(read.csv(expression_file, sep="\\t", header=TRUE, row.names=1)))\n  metadata <- as.matrix(read.csv(metadata_file, sep="\\t", header=TRUE, row.names=1))\n\n  print("Checking sample ordering...")\n  print(all.equal(colnames(expression_data), rownames(metadata)))\n\n  group <- interaction(metadata[,1])\n\n  mm <- model.matrix(~0 + group)\n\n  #print(head(expression_data))\n\n  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)\n\n  deseq_object <- DESeq(ddset)\n\n  deseq_results <- results(deseq_object, independentFiltering=FALSE)\n\n  deseq_results_df <-  as.data.frame(deseq_results)\n\n  # Save summary statistics of DEGs\n  if (data_type == "template") {\n    out_file = paste(local_dir, "DE_stats/DE_stats_template_data_", experiment_id,"_", run, ".txt", sep="")\n  } else if (data_type == "simulated") {\n    out_file = paste(local_dir, "DE_stats/DE_stats_simulated_data_", experiment_id,"_", run, ".txt", sep="")\n  }\n  write.table(deseq_results_df, file = out_file, row.names = T, sep = "\\t", quote = F)\n}\n\n\n# File created: "<local_dir>/DE_stats/DE_stats_template_data_SRP012656_real.txt"\nget_DE_stats_DESeq(metadata_filename,\n                   project_id, \n                   mapped_template_filename,\n                   "template",\n                   local_dir,\n                   "real_without_filtering")')


# In[40]:


template_nonrecount2_nofilter_DE_stats_filename = os.path.join(
    DE_stats_dir,
    "DE_stats_template_data_cis-par-KU1919_real_without_filtering.txt"
)

template_nonrecount2_nofilter_DE_stats = pd.read_csv(
    template_nonrecount2_nofilter_DE_stats_filename,
    sep="\t",
    index_col=0,
    header=0
)

template_nonrecount2_nofilter_DE_stats.loc["RBPJL"]


# In[41]:


template_nonrecount2_nofilter_DE_stats.loc["SSTR5.AS1"]


# **Takeaways:**
# * Case 1: genes with only simulated statistics is because template experiment has all 0 counts
# * Case 2: genes with fewer than 25 simulated experiments, some simulated experiments have all 0 counts for those genes
# * Case 3 (only found using non-recount2 template): genes with missing p-value in template experiment have NaN output from DESeq2, which indicates those genes that have been automatically filtered by DESeq2 as very likely to not have no significance: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
# 
# 
# **Proposed solution:**
# * Case 1 and 2: Remove genes with 0 counts across all samples. I will also add an option to additional remove genes that have mean counts < user specified threshold (for Jim's case). This should get rid of those rows with missing values in template statistic columns. These genes will be missing in the validation against Crow et. al. so I'll need to ignore these genes in the rank comparison. Any gene removed from the simulated experiments will get a lower `number of simulated experiments` reported so I will need to document this for the user.
# * Case 3: (option 1) I can change the parameter setting to turn off this autofiltering and it will perform all tests and report them. (option 2) Replace padj values = "NA" with "Filtered by DESeq2" and document that this means that DESeq2 pre-filtered these genes as likely not being signicant to help increase detection power. I'm not a super stats expert but I'm fine with option 2 if you are. DESeq2 documentation states: 
# ```
# filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic. Typically, this results in increased detection power at the same experiment-wide type I error, as measured in terms of the false discovery rate. 
# 
# For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from so high Poisson noise that any biological effect is drowned in the uncertainties from the read counting
# ```
