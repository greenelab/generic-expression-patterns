
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

# In[1]:


import os
import re
import pandas as pd
from ponyo import utils


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


# In[23]:


template_nonrecount2_DE_stats.loc["SSTR5.AS1"]


# **Takeaways:**
# * Case 1: genes with only simulated statistics is because template experiment has all 0 counts
# * Case 2: genes with fewer than 25 simulated experiments, some simulated experiments have all 0 counts for those genes
# * Case 3 (only found using non-recount2 template): genes with missing p-value in template experiment have NaN output from DESeq2, which indicates...
# 
# 
# **Proposed solution:**
# * Case 1: We can filter 0 count genes from template experiment to get rid of those rows with missing values in template statistic columns. These genes will be missing in the validation against Crow et. al. so I'll need to adjust for this.
# * Case 2: We can filter 0 count genes from simulated experiment but will still get a lower number of simulated experiments reported so this will need to be noted.
# * Case 3: 
