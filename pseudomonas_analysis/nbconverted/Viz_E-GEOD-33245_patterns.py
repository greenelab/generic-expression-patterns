
# coding: utf-8

# # Visualize E-GEOD-33245 patterns
# This notebook will examine patterns of generic and experiment-specific genes using E-GEOD-33245 as the template experiment
# 
# This experiment contains multiple comparisons/conditions:
# 
# * grp_1v2 compares WT vs *crc* mutants
# * grp_1v3 compares WT vs *cbrB* mutants
# * grp_1v4 compares WT vs *crcZ* mutant
# * grp_1v5 compares WT in LB vs WT in BSM

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import seaborn as sns
import numpy as np

from generic_expression_patterns_modules import process


# In[2]:


# Load data

# Summary data using abs value of test statistic
grp_1v2_file = "generic_gene_summary_E-GEOD-33245_1v2.tsv"
grp_1v3_file = "generic_gene_summary_E-GEOD-33245_1v3.tsv"
grp_1v4_file = "generic_gene_summary_E-GEOD-33245_1v4.tsv"
grp_1v5_file = "generic_gene_summary_E-GEOD-33245_1v5.tsv"

# Summary data using raw value of test statistic to get directionality
grp_1v2_raw_file = "generic_gene_summary_E-GEOD-33245_1v2_raw.tsv"
grp_1v3_raw_file = "generic_gene_summary_E-GEOD-33245_1v3_raw.tsv"
grp_1v4_raw_file = "generic_gene_summary_E-GEOD-33245_1v4_raw.tsv"
grp_1v5_raw_file = "generic_gene_summary_E-GEOD-33245_1v5_raw.tsv"


# In[3]:


# User parameters

# FDR adjusted p-value cutoff to use to define DEGs
pvalue_threshold = 0.05

# Get predicted generic DEGs using z-score cutoff
# Z-score cutoff was found by calculating the score
# whose invnorm(0.05/5549). Here we are using a p-value = 0.05
# with a Bonferroni correction for 5549 tests, which are
# the number of P. aeruginosa genes
zscore_threshold = 4.44


# ## Create dataframe to compare trends
# We are going to merge data across different conditions. For example, we will merge `grp_1v2` and `grp_1v3` to use for plotting later in this notebook. The Hogan lab can look at these tables to find *things of interest* as we start looking into how to use our computational predictions of generic and specific genes. 

# In[4]:


# Read data
grp_1v2 = pd.read_csv(grp_1v2_file, sep="\t", header=0, index_col=0)
grp_1v3 = pd.read_csv(grp_1v3_file, sep="\t", header=0, index_col=0)
grp_1v4 = pd.read_csv(grp_1v4_file, sep="\t", header=0, index_col=0)
grp_1v5 = pd.read_csv(grp_1v5_file, sep="\t", header=0, index_col=0)

grp_1v2_raw = pd.read_csv(grp_1v2_raw_file, sep="\t", header=0, index_col=0)
grp_1v3_raw = pd.read_csv(grp_1v3_raw_file, sep="\t", header=0, index_col=0)
grp_1v4_raw = pd.read_csv(grp_1v4_raw_file, sep="\t", header=0, index_col=0)
grp_1v5_raw = pd.read_csv(grp_1v5_raw_file, sep="\t", header=0, index_col=0)


# In[5]:


# Merge summary dfs using abs log2 FC and using raw values
merged_1v2s_df = process.merge_abs_raw_dfs(grp_1v2, grp_1v2_raw, '1v2')
merged_1v3s_df = process.merge_abs_raw_dfs(grp_1v3, grp_1v3_raw, '1v3')
merged_1v4s_df = process.merge_abs_raw_dfs(grp_1v4, grp_1v4_raw, '1v4')
merged_1v5s_df = process.merge_abs_raw_dfs(grp_1v5, grp_1v5_raw, '1v5')


# In[7]:


# Merge 1v2 and 1v3 summary dfs
merged_1v2_1v3_df = process.merge_two_conditions_df(merged_1v2s_df, merged_1v3s_df, '1v2', '1v3')
merged_1v2_1v3_df.head()


# In[6]:


# Merge 1v3 and 1v4 summaries
merged_1v3_1v4_df = process.merge_two_conditions_df(merged_1v3s_df, merged_1v4s_df, '1v3', '1v4')
merged_1v3_1v4_df.head()


# In[7]:


# Merge 1v2 and 1v4 summaries
merged_1v2_1v4_df = process.merge_two_conditions_df(merged_1v2s_df, merged_1v4s_df, '1v2', '1v4')
merged_1v2_1v4_df.head()


# In[8]:


# Save
merged_1v2_1v3_df.to_csv("merged_E-GEOD_1v2_1v3_directionality.tsv", sep="\t")
merged_1v2_1v4_df.to_csv("merged_E-GEOD_1v2_1v4_directionality.tsv", sep="\t")
merged_1v3_1v4_df.to_csv("merged_E-GEOD_1v3_1v4_directionality.tsv", sep="\t")


# ## Compare trends across different conditions
# 
# We want to compare across different conditions. For example, given:
# * Group 1v2: WT vs *crc* mutant
# * Group 1v3: WT vs *cbr* mutant
# 
# We can then compare the test statistic between these two groups above. We hope to see that,
# * Genes 1v3  > 1v2 are genes that change more in 1v3 than 1v2 and we guess are specific to 1v3 (high z-score)
# * Genes 1v3 < 1v2 are genes that change more in 1v2 than 1v3 and we guess are specific to 1v2 (high z-score)
# * Genes on the 1v3 = 1v2 line are those genes that change equally in both and we guess are generic genes (low z-score)

# ### 1v2 compared with 1v3

# In[9]:


process.plot_two_conditions(merged_1v2_1v3_df, "1v2", "1v3", "WT vs crc mutant", "WT vs cbrB mutant")


# ### 1v2 compared with 1v4

# In[10]:


process.plot_two_conditions(merged_1v2_1v4_df, "1v2", "1v4", "WT vs crc mutant", "WT vs crcZ mutant")


# ### 1v3 compared with 1v4

# In[11]:


process.plot_two_conditions(merged_1v3_1v4_df, "1v3", "1v4", "WT vs cbrB mutant", "WT vs crcZ mutant")


# **Takeaway:**
# Genes with high specificity score (i.e. genes with a high absolute value z-score) are located in the off x-y regions, as expected since these off diagonal regions represent those genes that are more differentially expressed in one condition versus the other. This shows some promise for using z-score to distinguish between generic and specific genes and we can start looking more into these trends.

# ## DEGs found using traditional criteria and using z-score
# 
# When performing DE analysis, this can return hundreds of genes using traditional criteria (FDR adjusted p-value < 0.05 and log2 fold change > 2). We want to see what genes are selected when we add z-score as an additional criteria to indicate genes that are specific to the pertubagen in question.

# ### 1v2

# In[10]:


(DEGs_1v2_df,
 degs_1v2_traditional,
 degs_1v2_specific,
 degs_1v2_generic,
 degs_1v2_intersect,
 degs_1v2_intersect_generic,
 degs_1v2_diff) = process.get_and_save_DEG_lists(merged_1v2s_df, '1v2', pvalue_threshold, zscore_threshold)


# In[11]:


process.plot_venn(degs_1v2_traditional, degs_1v2_specific, degs_1v2_generic)


# In[16]:


process.plot_volcanos(degs_1v2_intersect, degs_1v2_diff, merged_1v2s_df, "1v2", "WT vs crc mutant")


# ### 1v3

# In[15]:


(DEGs_1v3_df,
 degs_1v3_traditional,
 degs_1v3_specific,
 degs_1v3_generic,
 degs_1v3_intersect,
 degs_1v3_intersect_generic,
 degs_1v3_diff) = process.get_and_save_DEG_lists(merged_1v3s_df, '1v3', pvalue_threshold, zscore_threshold)


# In[16]:


process.plot_venn(degs_1v3_traditional, degs_1v3_specific, degs_1v3_generic)


# In[17]:


process.plot_volcanos(degs_1v3_intersect, degs_1v3_diff, merged_1v3s_df, "1v3", "WT vs cbrB mutant")


# ### 1v4

# In[18]:


(DEGs_1v4_df,
 degs_1v4_traditional,
 degs_1v4_specific,
 degs_1v4_generic,
 degs_1v4_intersect, 
 degs_1v4_intersect_generic,
 degs_1v4_diff) = process.get_and_save_DEG_lists(merged_1v4s_df, '1v4', pvalue_threshold, zscore_threshold)


# In[19]:


process.plot_venn(degs_1v4_traditional, degs_1v4_specific, degs_1v4_generic)


# In[20]:


process.plot_volcanos(degs_1v4_intersect, degs_1v4_diff, merged_1v4s_df, "1v4", "WT vs crcZ mutant")


# ### 1v5

# In[21]:


(DEGs_1v5_df,
 degs_1v5_traditional,
 degs_1v5_specific,
 degs_1v5_generic,
 degs_1v5_intersect,
 degs_1v5_intersect_generic,
 degs_1v5_diff) = process.get_and_save_DEG_lists(merged_1v5s_df, '1v5', pvalue_threshold, zscore_threshold)


# In[22]:


process.plot_venn(degs_1v5_traditional, degs_1v5_specific, degs_1v5_generic)


# In[23]:


process.plot_volcanos(degs_1v5_intersect, degs_1v5_diff, merged_1v5s_df, "1v5", "WT LB vs WT BSM")


# In[24]:


# Save DEGs to file to share with Hogan lab
degs_all_df = pd.DataFrame({'1v2 traditional + specific': pd.Series(degs_1v2_intersect),
                            '1v2 specific only': pd.Series(degs_1v2_diff),
                            '1v2 traditional + generic': pd.Series(degs_1v2_intersect_generic),
                            '1v3 traditional + specific': pd.Series(degs_1v3_intersect),
                            '1v3 specific only': pd.Series(degs_1v3_diff),
                            '1v3 traditional + generic': pd.Series(degs_1v3_intersect_generic),
                            '1v4 traditional + specific': pd.Series(degs_1v4_intersect),
                            '1v4 specific only': pd.Series(degs_1v4_diff),
                            '1v4 traditional + generic': pd.Series(degs_1v4_intersect_generic),
                            '1v5 traditional + specific': pd.Series(degs_1v5_intersect),
                            '1v5 specific only': pd.Series(degs_1v5_diff),
                            '1v5 traditional + generic': pd.Series(degs_1v5_intersect_generic)                            
                           })
degs_all_df.to_csv("DEGs_E_GEOD_33245.tsv", sep="\t")


# **Takeaway:**
# * Overall it looks like we find on the order of hundreds of DEGs using the traditional criteria (p-value and log2 fold change), but if we filter by high absolute z-score we get a reduced set of genes. In some cases the reduced set is of a similar size to the traditional set (i.e. 1v4) but other times the reduced set is signficantly smaller (i.e. 1v5).
# * Using traditional DE criteria you find a mix of generic and specific genes. And applying a z-score cutoff will supply researchers with a reasonable sized list of DEGs that are most relevant to the perturbagen in question that they can follow up with.
# * There are also some specific genes that are missed by the traditional DE criteria because they are just shy of being significant based on the FDR p-value cutoff. Perhaps these genes play an important niche role that would have otherwise been missed using traditional DE analysis.
