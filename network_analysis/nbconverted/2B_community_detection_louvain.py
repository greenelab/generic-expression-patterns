#!/usr/bin/env python
# coding: utf-8

# ## Run Louvain community detection experiments
# 
# To make sure that the results in `2A_community_detection_planted_partition.ipynb` aren't dependent on the choice of algorithm, this script runs the same analysis using the [Louvain community detection algorithm](https://en.wikipedia.org/wiki/Louvain_method). Here, we use the implementation in the [python-igraph package](https://igraph.org/python/doc/igraph.Graph-class.html#community_multilevel).
# 
# The Louvain algorithm optimizes [modularity](https://en.wikipedia.org/wiki/Louvain_method#Modularity_optimization), so in that sense it tends to find assortative (well-connected) communities similar to the PP model.

# In[1]:


import os
import random

import numpy as np
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
import seaborn as sns

from sample_nodes import (
    sample_degree_matched,
    sort_by_degree,
)

random.seed(1) # igraph uses Python RNG
np.random.seed(1)


# In[2]:


# relevant file paths
data_dir = './data'
edge_list = os.path.join(data_dir, 'edge_list_processed_unsigned.csv')

# map of Pa gene names to generic/not generic status, generated by Alex
generic_gene_map = os.path.join('..', 'pseudomonas_analysis', 'annot_df.tsv')

# script parameters
NUM_NODE_SAMPLES = 1000 # number of degree-matched node samples for permutation test
NUM_BINS = 100 # number of bins to divide nodes into, for sampling


# In[3]:


graph_df = pd.read_csv(edge_list)
graph_df.head()


# In[4]:


G = ig.Graph.TupleList(graph_df.values,
                       weights=True,
                       directed=False)


# In[5]:


# make sure vertex/edge properties exist
print(G.es['weight'][:5])


# In[6]:


annot_df = pd.read_csv(generic_gene_map, sep='\t', index_col=0)
annot_df.head()


# In[7]:


is_generic = [int(annot_df.loc[v['name'], 'label']) for v in G.vs]
G.vs['is_generic'] = is_generic


# In[8]:


# community detection using Louvain modularity optimization
partition = G.community_multilevel(weights=G.es['weight'], return_levels=False)
# plot?


# In[9]:


# get dataframe mapping Pa genes to communities
def graph_partition_to_df(G, partition):
    clusters = []
    for label, vl in enumerate(partition):
        clusters += [(G.vs['name'][v],
                      label,
                      G.degree(v),
                      G.vs['is_generic'][v]) for v in vl]
    return pd.DataFrame(clusters,
                        columns=['gene', 'label', 'degree', 'is_generic'])

labels_df = graph_partition_to_df(G, partition)
print(len(labels_df.label.unique()))
labels_df.sort_values(by='degree', ascending=False).head()


# In[10]:


# simultaneously sort nodes and degrees by degree, ascending
nodes, degrees, is_generic = sort_by_degree(labels_df.gene.values,
                                            labels_df.degree.values,
                                            labels_df.is_generic.values)

# sample a few times and add results to df
for it in range(NUM_NODE_SAMPLES):
    s_nodes, s_degrees, __ = sample_degree_matched(nodes, degrees, is_generic,
                                                   num_bins=NUM_BINS)
    sampled = [1 if gene in s_nodes else 0 for gene in labels_df.gene]
    labels_df['sampled_{}'.format(it)] = sampled

labels_df.sort_values(by='degree', ascending=False).iloc[:5, :5]


# In[11]:


generic_count_df = (
    labels_df.groupby('label').sum()
             .drop(columns=['degree'])
)
print(generic_count_df.shape)
generic_count_df.sort_values(by='is_generic', ascending=False).iloc[:5, :5]


# In[12]:


nonzero_counts_df = pd.DataFrame(
    [np.count_nonzero(generic_count_df, axis=0)],
    columns=generic_count_df.columns
)
nonzero_counts_df.iloc[:, :5]


# In[13]:


n_generic_groups = nonzero_counts_df.iloc[0, 0]
n_total_groups = len(labels_df.label.unique())

sns.set({'figure.figsize': (8, 5)})

sns.histplot(nonzero_counts_df.iloc[0, 1:].values,
             element='step',
             bins=np.arange(n_generic_groups, n_total_groups+1))
line1 = plt.gca().axvline(x=n_generic_groups,
                         linestyle='--', color='red')
line2 = plt.gca().axvline(x=n_total_groups,
                          linestyle='--', color='black')
plt.xlabel('Number of communities with selected genes')
plt.legend(handles=[mpatches.Patch(color=sns.color_palette()[0]), line1, line2],
           labels=['Random gene samples', 'Generic genes', 'Total communities'],
           loc='upper right')
plt.title('Number of communities, generic genes vs. {} degree-matched samples'.format(
    NUM_NODE_SAMPLES
))


# These results do not differ substantially from the planted partition results, suggesting that our conclusions are not sensitive to the choice of community detection algorithm.