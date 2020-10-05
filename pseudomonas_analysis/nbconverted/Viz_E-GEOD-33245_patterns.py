
# coding: utf-8

# # Visualize E-GEOD-33245 patterns
# This notebook will examine patterns of generic and experiment-specific genes using E-GEOD-33245 as the template experiment

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


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


# ## Create dataframe to compare trends
# We are going to merge data across different conditions. For example, we will merge `grp_1v2` and `grp_1v3` to use for plotting later in this notebook. Deb can look at these tables to find *things of interest* as we start looking into how to use our computational predictions of generic and specific genes. 
# 
# She would like a dataframe with 1v2 and 1v3, 1v3 and 1v4.

# In[3]:


# Read data
grp_1v2 = pd.read_csv(grp_1v2_file, sep="\t", header=0, index_col=0)
grp_1v3 = pd.read_csv(grp_1v3_file, sep="\t", header=0, index_col=0)
grp_1v4 = pd.read_csv(grp_1v4_file, sep="\t", header=0, index_col=0)
grp_1v5 = pd.read_csv(grp_1v5_file, sep="\t", header=0, index_col=0)

grp_1v2_raw = pd.read_csv(grp_1v2_raw_file, sep="\t", header=0, index_col=0)
grp_1v3_raw = pd.read_csv(grp_1v3_raw_file, sep="\t", header=0, index_col=0)
grp_1v4_raw = pd.read_csv(grp_1v4_raw_file, sep="\t", header=0, index_col=0)
grp_1v5_raw = pd.read_csv(grp_1v5_raw_file, sep="\t", header=0, index_col=0)


# In[ ]:


# FUNCTIONS
process.merge_one_condition()
process.merge_two_conditions_df(grp_1v2, grp_1v3, assoc dfs) --> returns merged df
process.plot_two_conditions(merged_df, grp_1v2, grp_1v3, xlabel, ylabel) --> returns plot
process.get_and save_DE_gene_lists(merged_one_condition) --> returns dict of lists
process.plot_volcanos(gene_lists, merged_one_conditions)
process.plot_venn(gene_lists,grp_1v2)


# In[4]:


# Merge 1v2 and 1v3 summaries
merged_1v2s_df = grp_1v2.merge(grp_1v2_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v2", "_grp_1v2_raw"])
merged_1v3s_df = grp_1v3.merge(grp_1v3_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v3", "_grp_1v3_raw"])
merged_1v2_1v3_all_df = merged_1v2s_df.merge(merged_1v3s_df, left_on='Gene ID', right_on="Gene ID")
merged_1v2_1v3_all_df.head()


# In[5]:


# Get specific columns requested by Deb
merged_1v2_1v3_all_df['max Z score'] = merged_1v2_1v3_all_df[['Z score_grp_1v2','Z score_grp_1v3']].max(axis=1)
merged_1v2_1v3_all_df['Gene ID Name'] = merged_1v2_1v3_all_df['Gene ID'] + " " + merged_1v2_1v3_all_df['Gene Name_grp_1v2'].fillna("")

merged_1v2_1v3_df = merged_1v2_1v3_all_df[['Gene ID',
                                           'Gene ID Name', 
                                           'Test statistic (Real)_grp_1v2',
                                           'Test statistic (Real)_grp_1v2_raw',
                                           'Adj P-value (Real)_grp_1v2',
                                           'Mean test statistic (simulated)_grp_1v2',
                                           'Std deviation (simulated)_grp_1v2',
                                           'Median adj p-value (simulated)_grp_1v2',
                                           'Test statistic (Real)_grp_1v3',
                                           'Test statistic (Real)_grp_1v3_raw',
                                           'Adj P-value (Real)_grp_1v3',
                                           'Mean test statistic (simulated)_grp_1v3',
                                           'Std deviation (simulated)_grp_1v3',
                                           'Median adj p-value (simulated)_grp_1v3',
                                           'Z score_grp_1v2',
                                           'Z score_grp_1v3',
                                           'max Z score'
                         ]]

merged_1v2_1v3_df.head()


# In[6]:


# Merge 1v3 and 1v4 summaries
merged_1v4s_df = grp_1v4.merge(grp_1v4_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v4", "_grp_1v4_raw"])
merged_1v3_1v4_all_df = merged_1v3s_df.merge(merged_1v4s_df, left_on='Gene ID', right_on="Gene ID")
merged_1v3_1v4_all_df.head()


# In[7]:


# Get specific columns requested by Deb
merged_1v3_1v4_all_df['max Z score'] = merged_1v3_1v4_all_df[['Z score_grp_1v3','Z score_grp_1v4']].max(axis=1)
merged_1v3_1v4_all_df['Gene ID Name'] = merged_1v3_1v4_all_df['Gene ID'] + " " + merged_1v3_1v4_all_df['Gene Name_grp_1v3'].fillna("")

merged_1v3_1v4_df = merged_1v3_1v4_all_df[['Gene ID',
                                           'Gene Name_grp_1v3', 
                                           'Test statistic (Real)_grp_1v3',
                                           'Test statistic (Real)_grp_1v3_raw',
                                           'Adj P-value (Real)_grp_1v3',
                                           'Mean test statistic (simulated)_grp_1v3',
                                           'Std deviation (simulated)_grp_1v3',
                                           'Median adj p-value (simulated)_grp_1v3',
                                           'Test statistic (Real)_grp_1v4',
                                           'Test statistic (Real)_grp_1v4_raw',
                                           'Adj P-value (Real)_grp_1v4',
                                           'Mean test statistic (simulated)_grp_1v4',
                                           'Std deviation (simulated)_grp_1v4',
                                           'Median adj p-value (simulated)_grp_1v4',
                                           'Z score_grp_1v3',
                                           'Z score_grp_1v4',
                                           'max Z score'
                         ]]

merged_1v3_1v4_df.head()


# In[8]:


# Merge 1v2 and 1v4 summaries
merged_1v2_1v4_all_df = merged_1v2s_df.merge(merged_1v4s_df, left_on='Gene ID', right_on="Gene ID")
merged_1v2_1v4_all_df.head()


# In[9]:


# Get specific columns requested by Deb
merged_1v2_1v4_all_df['max Z score'] = merged_1v2_1v4_all_df[['Z score_grp_1v2','Z score_grp_1v4']].max(axis=1)
merged_1v2_1v4_all_df['Gene ID Name'] = merged_1v2_1v4_all_df['Gene ID'] + " " + merged_1v2_1v4_all_df['Gene Name_grp_1v2'].fillna("")

merged_1v2_1v4_df = merged_1v2_1v4_all_df[['Gene ID',
                                           'Gene Name_grp_1v2', 
                                           'Test statistic (Real)_grp_1v2',
                                           'Test statistic (Real)_grp_1v2_raw',
                                           'Adj P-value (Real)_grp_1v2',
                                           'Mean test statistic (simulated)_grp_1v2',
                                           'Std deviation (simulated)_grp_1v2',
                                           'Median adj p-value (simulated)_grp_1v2',
                                           'Test statistic (Real)_grp_1v4',
                                           'Test statistic (Real)_grp_1v4_raw',
                                           'Adj P-value (Real)_grp_1v4',
                                           'Mean test statistic (simulated)_grp_1v4',
                                           'Std deviation (simulated)_grp_1v4',
                                           'Median adj p-value (simulated)_grp_1v4',
                                           'Z score_grp_1v2',
                                           'Z score_grp_1v4',
                                           'max Z score'
                         ]]

merged_1v2_1v4_df.head()


# In[10]:


# Save
merged_1v2_1v3_df.to_csv("merged_E-GEOD_1v2_1v3_directionality.tsv", sep="\t")
merged_1v2_1v4_df.to_csv("merged_E-GEOD_1v2_1v4_directionality.tsv", sep="\t")
merged_1v3_1v4_df.to_csv("merged_E-GEOD_1v3_1v4_directionality.tsv", sep="\t")


# ## Compare trends across different conditions
# 
# We want to compare across different conditions. For example, given:
# * Group 1v2: WT vs crc mutant
# * Group 1v3: WT vs cbr mutant
# 
# We can then compare the test statistic between these two groups above. We hope to see that,
# * Genes 1v3  > 1v2 are genes that change more in 1v3 than 1v2 and we guess are specific to 1v3 (high z-score)
# * Genes 1v3 < 1v2 are genes that change more in 1v2 than 1v3 and we guess are specific to 1v2 (high z-score)
# * Genes on the 1v3 = 1v2 line are those genes that change equally in both and we guess are generic genes (low z-score)

# ### 1v2 compared with 1v3

# In[11]:


fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)

fig_1v2_abs = sns.scatterplot(data=merged_1v2_1v3_df,
                              x="Test statistic (Real)_grp_1v2",
                              y="Test statistic (Real)_grp_1v3",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[0],
                              palette=cmap)
fig_1v2_abs.plot([0,4],[0,4],"--k")

fig_1v2_raw = sns.scatterplot(data=merged_1v2_1v3_df,
                              x="Test statistic (Real)_grp_1v2_raw",
                              y="Test statistic (Real)_grp_1v3_raw",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[1],
                              palette=cmap)
fig_1v2_raw.plot([-4,4],[-4,4],"--k")

# Add labels
fig.suptitle('crc mutant vs cbrB mutant', fontsize=16)
fig.text(0.5, 0.04, 'WT vs crc mutant', ha='center', va='center')
fig.text(0.06, 0.5, 'WT vs cbrB mutant', ha='center', va='center', rotation='vertical')
axes[0].set_title('using abs(log2FC)')
axes[1].set_title('using log2FC')
axes[0].set_xlabel('')
axes[1].set_xlabel('')
axes[0].set_ylabel('')
axes[1].set_ylabel('')
print(fig)

# ADD NEWLINE TO TITLE
# MAKE LABELS LARGER
# MOVE LEGEND?


# ### 1v2 compared with 1v4

# In[12]:


fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)

fig_1v2_abs = sns.scatterplot(data=merged_1v2_1v4_df,
                              x="Test statistic (Real)_grp_1v2",
                              y="Test statistic (Real)_grp_1v4",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[0],
                              palette=cmap)
fig_1v2_abs.plot([0,4],[0,4],"--k")

fig_1v2_raw = sns.scatterplot(data=merged_1v2_1v4_df,
                              x="Test statistic (Real)_grp_1v2_raw",
                              y="Test statistic (Real)_grp_1v4_raw",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[1],
                              palette=cmap)
fig_1v2_raw.plot([-4,4],[-4,4],"--k")

# Add labels
fig.suptitle('crc mutant vs crcZ mutant', fontsize=16)
fig.text(0.5, 0.04, 'WT vs crc mutant', ha='center', va='center')
fig.text(0.06, 0.5, 'WT vs crcZ mutant', ha='center', va='center', rotation='vertical')
axes[0].set_title('using abs(log2FC)')
axes[1].set_title('using log2FC')
axes[0].set_xlabel('')
axes[1].set_xlabel('')
axes[0].set_ylabel('')
axes[1].set_ylabel('')
print(fig)


# ### 1v3 compared with 1v4

# In[13]:


fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)

fig_1v2_abs = sns.scatterplot(data=merged_1v3_1v4_df,
                              x="Test statistic (Real)_grp_1v3",
                              y="Test statistic (Real)_grp_1v4",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[0],
                              palette=cmap)
fig_1v2_abs.plot([0,4],[0,4],"--k")

fig_1v2_raw = sns.scatterplot(data=merged_1v3_1v4_df,
                              x="Test statistic (Real)_grp_1v3_raw",
                              y="Test statistic (Real)_grp_1v4_raw",
                              hue="max Z score",
                              size="max Z score",
                              linewidth=0,
                              alpha=0.7,
                              ax=axes[1],
                              palette=cmap)
fig_1v2_raw.plot([-4,4],[-4,4],"--k")

# Add labels
fig.suptitle('cbrB mutant vs crcZ mutant', fontsize=16)
fig.text(0.5, 0.04, 'WT vs cbrB mutant', ha='center', va='center')
fig.text(0.06, 0.5, 'WT vs crcZ mutant', ha='center', va='center', rotation='vertical')
axes[0].set_title('using abs(log2FC)')
axes[1].set_title('using log2FC')
axes[0].set_xlabel('')
axes[1].set_xlabel('')
axes[0].set_ylabel('')
axes[1].set_ylabel('')
print(fig)


# **Takeaway:**
# * A few specific genes (with very high z-score) in the off x-y regions as expected. This shows some promise for using z-score to distinguish between generic and specific genes and we can start looking more into these trends.
# 
# ADD TEXT HERE

# ## DEGs found using traditional criteria and using z-score
# 
# When performing DE analysis, this can return hundreds of genes using traditional criteria (FDR adjusted p-value < 0.05 and log2 fold change > 2). We want to see what genes are selected when we add z-score as an additional criteria to indicate genes that are specific to the pertubagen in question.

# ### 1v2

# In[14]:


# Get DEGs using traditional criteria
degs_1v2_traditional = list((merged_1v2s_df[(merged_1v2s_df['Test statistic (Real)_grp_1v2']>1)
                                           & (merged_1v2s_df['Adj P-value (Real)_grp_1v2']<0.05)]
                             .set_index('Gene ID')
                             .index)
                           )
print(f'No. of DEGs using traditional criteria: {len(degs_1v2_traditional)}')

# Get predicted specific DEGs using z-score cutoff
degs_1v2_specific = list((merged_1v2s_df[(merged_1v2s_df['Test statistic (Real)_grp_1v2']>1)
                                           & (merged_1v2s_df['Z score_grp_1v2']>10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of specific DEGs using z-score: {len(degs_1v2_specific)}')

# Get predicted generic DEGs using z-score cutoff
degs_1v2_generic = list((merged_1v2s_df[(merged_1v2s_df['Test statistic (Real)_grp_1v2']>1)
                                           & (merged_1v2s_df['Z score_grp_1v2']<10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of generic DEGs using z-score: {len(degs_1v2_generic)}')

# Get intersection of DEGs using traditional and z-score criteria
degs_1v2_intersect = list(set(degs_1v2_traditional).intersection(degs_1v2_specific))
print(f'No. of traditional DEGs that pass z-score criteria: {len(degs_1v2_intersect)}')

# Get specific DEGs that were NOT found using traditional criteria
degs_1v2_diff = list(set(degs_1v2_specific).difference(degs_1v2_intersect))
print(f'No. of specific DEGs that were not found by traditional criteria: {len(degs_1v2_diff)}')


# In[15]:


# Save list of genes that interesect and those that do not
merged_1v2s_df['Gene ID Name'] = merged_1v2s_df['Gene ID'] + " " + merged_1v2s_df['Gene Name_grp_1v2'].fillna("")

# Set `Gene ID` as index
merged_1v2s_df.set_index('Gene ID', inplace=True)

gene_id_names_1v2_intersect = merged_1v2s_df.loc[degs_1v2_intersect, 'Gene ID Name']
gene_id_names_1v2_diff = merged_1v2s_df.loc[degs_1v2_diff, 'Gene ID Name']
gene_id_names_1v2_generic = merged_1v2s_df.loc[degs_1v2_generic, 'Gene ID Name']

gene_lists_1v2_df = pd.DataFrame({'Traditional + specific DEGs': gene_id_names_1v2_intersect,
                                  'Specific only DEGs': gene_id_names_1v2_diff,
                                  'Generic DEGs': gene_id_names_1v2_generic
                                 }
                                )


# In[16]:


fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))

# Add columns for plotting
merged_1v2s_df['FDR adjuted p-value plot'] = -np.log10(merged_1v2s_df['Adj P-value (Real)_grp_1v2'])
merged_1v2s_df['gene group'] = 'none'
merged_1v2s_df.loc[degs_1v2_intersect, 'gene group'] = 'traditional + specific DEGs'
merged_1v2s_df.loc[degs_1v2_diff,'gene group'] = 'specific only DEGs'

# Plot: log2FC vs p-value (traditional criteria)
sns.scatterplot(data=merged_1v2s_df,
                x="Test statistic (Real)_grp_1v2_raw",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[0])

# Plot: log2FC vs z-score
sns.scatterplot(data=merged_1v2s_df,
                x="Test statistic (Real)_grp_1v2_raw",
                y="Z score_grp_1v2",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[1])

# Plot: z-score vs p-value
sns.scatterplot(data=merged_1v2s_df,
                x="Z score_grp_1v2",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[2])

# Add labels
fig.suptitle('WT vs crc mutant', fontsize=16)
axes[0].set_xlabel('log2 Fold Change')
axes[1].set_xlabel('log2 Fold Change')
axes[2].set_xlabel('Z-score')
axes[0].set_ylabel('FDR adjusted p-value')
axes[1].set_ylabel('Z-score (specificity)')
axes[2].set_ylabel('FDR adjusted p-value')
axes[0].set_title('log2FC vs p-value')
axes[1].set_title('log2FC vs z-score')
print(fig)

# ADD NEWLINE TO TITLE
# ADD THRESHOLD
# ADD P-VALUE VS Z-SCORE


# ### 1v3

# In[17]:


# Get DEGs using traditional criteria
degs_1v3_traditional = list((merged_1v3s_df[(merged_1v3s_df['Test statistic (Real)_grp_1v3']>1)
                                           & (merged_1v3s_df['Adj P-value (Real)_grp_1v3']<0.05)]
                             .set_index('Gene ID')
                             .index)
                           )
print(f'No. of DEGs using traditional criteria: {len(degs_1v3_traditional)}')

# Get predicted specific DEGs using z-score cutoff
degs_1v3_specific = list((merged_1v3s_df[(merged_1v3s_df['Test statistic (Real)_grp_1v3']>1)
                                           & (merged_1v3s_df['Z score_grp_1v3']>10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of specific DEGs using z-score: {len(degs_1v2_specific)}')

# Get predicted generic DEGs using z-score cutoff
degs_1v3_generic = list((merged_1v3s_df[(merged_1v3s_df['Test statistic (Real)_grp_1v3']>1)
                                           & (merged_1v3s_df['Z score_grp_1v3']<10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of generic DEGs using z-score: {len(degs_1v3_generic)}')

# Get intersection of DEGs using traditional and z-score criteria
degs_1v3_intersect = list(set(degs_1v3_traditional).intersection(degs_1v3_specific))
print(f'No. of traditional DEGs that pass z-score criteria: {len(degs_1v3_intersect)}')

# Get specific DEGs that were NOT found using traditional criteria
degs_1v3_diff = list(set(degs_1v3_specific).difference(degs_1v3_intersect))
print(f'No. of specific DEGs that were not found by traditional criteria: {len(degs_1v3_diff)}')


# In[18]:


# Save list of genes that interesect and those that do not
merged_1v3s_df['Gene ID Name'] = merged_1v3s_df['Gene ID'] + " " + merged_1v3s_df['Gene Name_grp_1v3'].fillna("")

# Set `Gene ID` as index
merged_1v3s_df.set_index('Gene ID', inplace=True)

gene_id_names_1v3_intersect = merged_1v2s_df.loc[degs_1v3_intersect, 'Gene ID Name']
gene_id_names_1v3_diff = merged_1v2s_df.loc[degs_1v3_diff, 'Gene ID Name']
gene_id_names_1v3_generic = merged_1v2s_df.loc[degs_1v3_generic, 'Gene ID Name']

gene_lists_1v3_df = pd.DataFrame({'Traditional + specific DEGs': gene_id_names_1v3_intersect,
                                  'Specific only DEGs': gene_id_names_1v3_diff,
                                  'Generic DEGs': gene_id_names_1v3_generic
                                 }
                                )


# In[19]:


fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))

# Add columns for plotting
merged_1v3s_df['FDR adjuted p-value plot'] = -np.log10(merged_1v3s_df['Adj P-value (Real)_grp_1v3'])
merged_1v3s_df['gene group'] = 'none'
merged_1v3s_df.loc[degs_1v3_intersect, 'gene group'] = 'traditional + specific DEGs'
merged_1v3s_df.loc[degs_1v3_diff,'gene group'] = 'specific only DEGs'

# Plot: log2FC vs p-value (traditional criteria)
sns.scatterplot(data=merged_1v3s_df,
                x="Test statistic (Real)_grp_1v3_raw",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[0])

# Plot: log2FC vs z-score
sns.scatterplot(data=merged_1v3s_df,
                x="Test statistic (Real)_grp_1v3_raw",
                y="Z score_grp_1v3",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[1])

# Plot: z-score vs p-value
sns.scatterplot(data=merged_1v3s_df,
                x="Z score_grp_1v3",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[2])

# Add labels
fig.suptitle('WT vs cbrB mutant', fontsize=16)
axes[0].set_xlabel('log2 Fold Change')
axes[1].set_xlabel('log2 Fold Change')
axes[2].set_xlabel('Z-score')
axes[0].set_ylabel('FDR adjusted p-value')
axes[1].set_ylabel('Z-score (specificity)')
axes[2].set_ylabel('FDR adjusted p-value')
axes[0].set_title('log2FC vs p-value')
axes[1].set_title('log2FC vs z-score')
print(fig)


# ### 1v4

# In[20]:


# Get DEGs using traditional criteria
degs_1v4_traditional = list((merged_1v4s_df[(merged_1v4s_df['Test statistic (Real)_grp_1v4']>1)
                                           & (merged_1v4s_df['Adj P-value (Real)_grp_1v4']<0.05)]
                             .set_index('Gene ID')
                             .index)
                           )
print(f'No. of DEGs using traditional criteria: {len(degs_1v4_traditional)}')

# Get predicted specific DEGs using z-score cutoff
degs_1v4_specific = list((merged_1v4s_df[(merged_1v4s_df['Test statistic (Real)_grp_1v4']>1)
                                           & (merged_1v4s_df['Z score_grp_1v4']>10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of specific DEGs using z-score: {len(degs_1v4_specific)}')

# Get predicted generic DEGs using z-score cutoff
degs_1v4_generic = list((merged_1v4s_df[(merged_1v4s_df['Test statistic (Real)_grp_1v4']>1)
                                           & (merged_1v4s_df['Z score_grp_1v4']<10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of generic DEGs using z-score: {len(degs_1v4_generic)}')

# Get intersection of DEGs using traditional and z-score criteria
degs_1v4_intersect = list(set(degs_1v4_traditional).intersection(degs_1v4_specific))
print(f'No. of traditional DEGs that pass z-score criteria: {len(degs_1v4_intersect)}')

# Get specific DEGs that were NOT found using traditional criteria
degs_1v4_diff = list(set(degs_1v4_specific).difference(degs_1v4_intersect))
print(f'No. of specific DEGs that were not found by traditional criteria: {len(degs_1v4_diff)}')


# In[21]:


# Save list of genes that interesect and those that do not
merged_1v4s_df['Gene ID Name'] = merged_1v4s_df['Gene ID'] + " " + merged_1v4s_df['Gene Name_grp_1v4'].fillna("")

# Set `Gene ID` as index
merged_1v4s_df.set_index('Gene ID', inplace=True)

gene_id_names_1v4_intersect = merged_1v4s_df.loc[degs_1v4_intersect, 'Gene ID Name']
gene_id_names_1v4_diff = merged_1v4s_df.loc[degs_1v4_diff, 'Gene ID Name']
gene_id_names_1v4_generic = merged_1v4s_df.loc[degs_1v4_generic, 'Gene ID Name']

gene_lists_1v4_df = pd.DataFrame({'Traditional + specific DEGs': gene_id_names_1v4_intersect,
                                  'Specific only DEGs': gene_id_names_1v4_diff,
                                  'Generic DEGs': gene_id_names_1v4_generic
                                 }
                                )


# In[22]:


fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))

# Add columns for plotting
merged_1v4s_df['FDR adjuted p-value plot'] = -np.log10(merged_1v4s_df['Adj P-value (Real)_grp_1v4'])
merged_1v4s_df['gene group'] = 'none'
merged_1v4s_df.loc[degs_1v4_intersect, 'gene group'] = 'traditional + specific DEGs'
merged_1v4s_df.loc[degs_1v4_diff,'gene group'] = 'specific only DEGs'

# Plot: log2FC vs p-value (traditional criteria)
sns.scatterplot(data=merged_1v4s_df,
                x="Test statistic (Real)_grp_1v4_raw",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[0])

# Plot: log2FC vs z-score
sns.scatterplot(data=merged_1v4s_df,
                x="Test statistic (Real)_grp_1v4_raw",
                y="Z score_grp_1v4",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[1])

# Plot: z-score vs p-value
sns.scatterplot(data=merged_1v4s_df,
                x="Z score_grp_1v4",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.7,
                ax=axes[2])

# Add labels
fig.suptitle('WT vs crcZ mutant', fontsize=16)
axes[0].set_xlabel('log2 Fold Change')
axes[1].set_xlabel('log2 Fold Change')
axes[2].set_xlabel('Z-score')
axes[0].set_ylabel('FDR adjusted p-value')
axes[1].set_ylabel('Z-score (specificity)')
axes[2].set_ylabel('FDR adjusted p-value')
axes[0].set_title('log2FC vs p-value')
axes[1].set_title('log2FC vs z-score')
print(fig)


# ### 1v5

# In[23]:


# Create merged df using grp_1v5_raw and grp-1v5
merged_1v5s_df = grp_1v5.merge(grp_1v5_raw, left_on='Gene ID', right_on="Gene ID", suffixes=["_grp_1v5", "_grp_1v5_raw"])

# Get DEGs using traditional criteria
degs_1v5_traditional = list((merged_1v5s_df[(merged_1v5s_df['Test statistic (Real)_grp_1v5']>1)
                                           & (merged_1v5s_df['Adj P-value (Real)_grp_1v5']<0.05)]
                             .set_index('Gene ID')
                             .index)
                           )
print(f'No. of DEGs using traditional criteria: {len(degs_1v5_traditional)}')

# Get predicted specific DEGs using z-score cutoff
degs_1v5_specific = list((merged_1v5s_df[(merged_1v5s_df['Test statistic (Real)_grp_1v5']>1)
                                           & (merged_1v5s_df['Z score_grp_1v5']>10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of specific DEGs using z-score: {len(degs_1v5_specific)}')

# Get predicted generic DEGs using z-score cutoff
degs_1v5_generic = list((merged_1v5s_df[(merged_1v5s_df['Test statistic (Real)_grp_1v5']>1)
                                           & (merged_1v5s_df['Z score_grp_1v5']<10)]
                          .set_index('Gene ID')
                          .index)
                           )
print(f'No. of generic DEGs using z-score: {len(degs_1v5_generic)}')

# Get intersection of DEGs using traditional and z-score criteria
degs_1v5_intersect = list(set(degs_1v5_traditional).intersection(degs_1v5_specific))
print(f'No. of traditional DEGs that pass z-score criteria: {len(degs_1v5_intersect)}')

# Get specific DEGs that were NOT found using traditional criteria
degs_1v5_diff = list(set(degs_1v5_specific).difference(degs_1v5_intersect))
print(f'No. of specific DEGs that were not found by traditional criteria: {len(degs_1v5_diff)}')


# In[24]:


# Save list of genes that interesect and those that do not
merged_1v5s_df['Gene ID Name'] = merged_1v5s_df['Gene ID'] + " " + merged_1v5s_df['Gene Name_grp_1v5'].fillna("")

# Set `Gene ID` as index
merged_1v5s_df.set_index('Gene ID', inplace=True)

gene_id_names_1v5_intersect = merged_1v5s_df.loc[degs_1v5_intersect, 'Gene ID Name']
gene_id_names_1v5_diff = merged_1v5s_df.loc[degs_1v5_diff, 'Gene ID Name']
gene_id_names_1v5_generic = merged_1v5s_df.loc[degs_1v5_generic, 'Gene ID Name']

gene_lists_1v5_df = pd.DataFrame({'Traditional + specific DEGs': gene_id_names_1v5_intersect,
                                  'Specific only DEGs': gene_id_names_1v5_diff,
                                  'Generic DEGs': gene_id_names_1v5_generic
                                 }
                                )


# In[25]:


fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))

# Add columns for plotting
merged_1v5s_df['FDR adjuted p-value plot'] = -np.log10(merged_1v5s_df['Adj P-value (Real)_grp_1v5'])
merged_1v5s_df['gene group'] = 'none'
merged_1v5s_df.loc[degs_1v2_intersect, 'gene group'] = 'traditional + specific DEGs'
merged_1v5s_df.loc[degs_1v2_diff,'gene group'] = 'specific only DEGs'

# Plot: log2FC vs p-value (traditional criteria)
sns.scatterplot(data=merged_1v5s_df,
                x="Test statistic (Real)_grp_1v5_raw",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[0])

# Plot: log2FC vs z-score
sns.scatterplot(data=merged_1v5s_df,
                x="Test statistic (Real)_grp_1v5_raw",
                y="Z score_grp_1v5",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[1])

# Plot: z-score vs p-value
sns.scatterplot(data=merged_1v5s_df,
                x="Z score_grp_1v5",
                y="FDR adjuted p-value plot",
                hue="gene group",
                linewidth=0,
                alpha=0.5,
                ax=axes[2])

# Add labels
fig.suptitle('WT in LB vs WT in BSM', fontsize=16)
axes[0].set_xlabel('log2 Fold Change')
axes[1].set_xlabel('log2 Fold Change')
axes[2].set_xlabel('Z-score')
axes[0].set_ylabel('FDR adjusted p-value')
axes[1].set_ylabel('Z-score (specificity)')
axes[2].set_ylabel('FDR adjusted p-value')
axes[0].set_title('log2FC vs p-value')
axes[1].set_title('log2FC vs z-score')
print(fig)


# **Takeaway:**
# * Overall it looks like we find on the order of hundreds of DEGs using the traditional criteria (p-value and log2 fold change), but if we filter by z-score we get a reduced set of genes
# * We predict that this z-score cutoff will supply researchers with a reasonable sized list of DEGs to follow-up with and that these are the genes that are most relevant to the perturbagen in question.
