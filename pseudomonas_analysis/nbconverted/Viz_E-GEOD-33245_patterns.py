
# coding: utf-8

# # Visualize E-GEOD-33245 patterns
# This notebook will examine patterns of generic and experiment-specific genes using E-GEOD-33245 as the template experiment

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


# Load data
# Summary data using abs value of test statistic
grp_1v2_file = "generic_gene_summary_E-GEOD-33245_1v2.tsv"
grp_1v3_file = "generic_gene_summary_E-GEOD-33245_1v3.tsv"
grp_1v4_file = "generic_gene_summary_E-GEOD-33245_1v4.tsv"

# Summary data using raw value of test statistic to get directionality
grp_1v2_raw_file = "generic_gene_summary_E-GEOD-33245_1v2_raw.tsv"
grp_1v3_raw_file = "generic_gene_summary_E-GEOD-33245_1v3_raw.tsv"
grp_1v4_raw_file = "generic_gene_summary_E-GEOD-33245_1v4_raw.tsv"


# ## Compare 1v2 and 1v3 grouping
# Within experiment E-GEOD-33245 we can perform 2 comparisons:
# * Group 1v2: WT vs crc mutant
# * Group 1v3: WT vs cbr mutant
# 
# We can then compare the test statistic between these two groups above. We hope to see that,
# * Genes 1v3  > 1v2 are genes that change more in 1v3 than 1v2 and we guess are specific to 1v3 (high z-score)
# * Genes 1v3 < 1v2 are genes that change more in 1v2 than 1v3 and we guess are specific to 1v2 (high z-score)
# * Genes on the 1v3 = 1v2 line are those genes that change equally in both and we guess are generic genes (low z-score)

# In[3]:


# Read data
grp_1v2 = pd.read_csv(grp_1v2_file, sep="\t", header=0, index_col=0)
grp_1v3 = pd.read_csv(grp_1v3_file, sep="\t", header=0, index_col=0)


# In[4]:


grp_1v2.head()


# In[5]:


grp_1v3.head()


# In[6]:


merged_df = grp_1v2.merge(grp_1v3, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v2", "_grp_1v3"])


# In[7]:


merged_df.head()


# In[8]:


# Add max(z-score of 1v2, z-score 1v3) to color by
merged_df['mean Z score'] = abs(merged_df[['Z score_grp_1v2','Z score_grp_1v3']].mean(axis=1))
merged_df['max Z score'] = abs(merged_df[['Z score_grp_1v2','Z score_grp_1v3']].max(axis=1))
merged_df.head()


# In[9]:


sns.distplot(merged_df['max Z score'], kde=False)


# In[10]:


cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)
sns.scatterplot(data=merged_df,
                #data=merged_df[merged_df["max Z score"]<0.5],
                x="Test statistic (Real)_grp_1v2",
                y="Test statistic (Real)_grp_1v3",
                hue="max Z score",
                size="max Z score",
                linewidth=0,
                alpha=0.7,
                palette=cmap)
plt.plot([0,4],[0,4],"--k")


# In[11]:


merged_df[merged_df['max Z score']>20]


# In[12]:


merged_df.to_csv("merged_E-GEOD_1v2_1v3.tsv", sep="\t")


# **Takeaway:**
# * A few specific genes (with very high z-score) in the off x-y regions as expected. This shows some promise for using z-score to distinguish between generic and specific genes and we can start looking more into these trends.

# ## Create dataframe to compare trends
# Here we create dataframes that Deb can look at to find *things of interest* as we start looking into how to use our computational predictions of generic and specific genes.
# 
# She would like a dataframe with 1v2 and 1v3, 1v3 and 1v4.

# In[13]:


# Read data
grp_1v2 = pd.read_csv(grp_1v2_file, sep="\t", header=0, index_col=0)
grp_1v3 = pd.read_csv(grp_1v3_file, sep="\t", header=0, index_col=0)
grp_1v4 = pd.read_csv(grp_1v4_file, sep="\t", header=0, index_col=0)

grp_1v2_raw = pd.read_csv(grp_1v2_raw_file, sep="\t", header=0, index_col=0)
grp_1v3_raw = pd.read_csv(grp_1v3_raw_file, sep="\t", header=0, index_col=0)
grp_1v4_raw = pd.read_csv(grp_1v4_raw_file, sep="\t", header=0, index_col=0)


# In[14]:


# Merge 1v2 and 1v3 summaries
merged_1v2s_df = grp_1v2.merge(grp_1v2_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v2", "_grp_1v2_raw"])
merged_1v3s_df = grp_1v3.merge(grp_1v3_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v3", "_grp_1v3_raw"])
merged_1v2_1v3_all_df = merged_1v2s_df.merge(merged_1v3s_df, left_on='Gene ID', right_on="Gene ID")
merged_1v2_1v3_all_df.head()


# In[15]:


# Get specific columns requested by Deb
merged_1v2_1v3_all_df['max Z score'] = merged_1v2_1v3_all_df[['Z score_grp_1v2','Z score_grp_1v3']].max(axis=1)

merged_1v2_1v3_df = merged_1v2_1v3_all_df[['Gene ID',
                                           'Gene Name_grp_1v2', 
                                           'Test statistic (Real)_grp_1v2',
                                           'Test statistic (Real)_grp_1v2_raw',
                                           'Test statistic (Real)_grp_1v3',
                                           'Test statistic (Real)_grp_1v3_raw',
                                           'Z score_grp_1v2',
                                           'Z score_grp_1v3',
                                           'max Z score'
                         ]]

merged_1v2_1v3_df.head()


# In[16]:


# Merge 1v3 and 1v4 summaries
merged_1v4s_df = grp_1v4.merge(grp_1v4_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v4", "_grp_1v4_raw"])
merged_1v3_1v4_all_df = merged_1v3s_df.merge(merged_1v4s_df, left_on='Gene ID', right_on="Gene ID")
merged_1v3_1v4_all_df.head()


# In[17]:


# Get specific columns requested by Deb
merged_1v3_1v4_all_df['max Z score'] = merged_1v3_1v4_all_df[['Z score_grp_1v3','Z score_grp_1v4']].max(axis=1)
merged_1v3_1v4_df = merged_1v3_1v4_all_df[['Gene ID',
                                           'Gene Name_grp_1v3', 
                                           'Test statistic (Real)_grp_1v3',
                                           'Test statistic (Real)_grp_1v3_raw',
                                           'Test statistic (Real)_grp_1v4',
                                           'Test statistic (Real)_grp_1v4_raw',
                                           'Z score_grp_1v3',
                                           'Z score_grp_1v4',
                                           'max Z score'
                         ]]

merged_1v3_1v4_df.head()


# In[18]:


# Save
merged_1v2_1v3_df.to_csv("merged_E-GEOD_1v2_1v3_directionality.tsv", sep="\t")
merged_1v3_1v4_df.to_csv("merged_E-GEOD_1v3_1v4_directionality.tsv", sep="\t")

