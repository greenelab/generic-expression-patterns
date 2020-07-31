
# coding: utf-8

# # Visualize E-GEOD-33245 patterns
# This notebook will examine patterns of generic and experiment-specific genes using E-GEOD-33245 as the template experiment

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd


# In[2]:


# Load data
grp_1v2_file = "generic_gene_summary_E-GEOD-33245_1v2.tsv"
grp_1v3_file = "generic_gene_summary_E-GEOD-33245_1v3.tsv"


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


# In[21]:


# Add max(z-score of 1v2, z-score 1v3) to color by
merged_df['max Z score'] = abs(merged_df[['Z score_grp_1v2','Z score_grp_1v3']].max(axis=1))
merged_df.head()


# In[34]:


import seaborn as sns
import matplotlib.pyplot as plt
sns.scatterplot(data=merged_df,
                #data=merged_df[merged_df["max Z score"]<0.5],
                x="Mean test statistic (simulated)_grp_1v2",
                y="Mean test statistic (simulated)_grp_1v3",
                hue="max Z score",
                size="max Z score",
                alpha=1)
plt.plot([0,0.3],[0,0.3])


# In[13]:


merged_df[merged_df['max Z score']>20]


# In[10]:


#merged_df.to_csv("merged_E-GEOD_1v2_1v3.tsv", sep="\t")

