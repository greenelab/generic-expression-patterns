#!/usr/bin/env python
# coding: utf-8

# # Visualize recount2 compendium
# This notebook does a quick visualization of the recount2 compendium (subset of recount2) being used in training the VAE model just to get a sense for the amount of variation in the data and to roughly determine if the training and validation set have similar distributions.
# 
# This will help to interpret the loss curves.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd
import numpy as np
import random
import umap

from plotnine import (ggplot,
                      labs,  
                      geom_line, 
                      geom_point,
                      geom_errorbar,
                      aes, 
                      ggsave, 
                      theme_bw,
                      theme,
                      xlim,
                      ylim,
                      facet_wrap,
                      scale_color_manual,
                      guides, 
                      guide_legend,
                      element_blank,
                      element_text,
                      element_rect,
                      element_line,
                      coords)

from ponyo import utils

from numpy.random import seed
random_state = 123
seed(random_state)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

config_file = os.path.abspath(os.path.join(base_dir,
                                           "config_human.tsv"))
params = utils.read_config(config_file)


# In[3]:


# Load params
normalized_data_file = params['normalized_compendium_data_file']
validation_frac = params['validation_frac']


# In[4]:


# Read data
normalized_compendium = pd.read_table(
    normalized_data_file,
    header=0,
    sep='\t',
    index_col=0)

print(normalized_compendium.shape)
normalized_compendium.head()


# In[5]:


# Get validation and training set
# random_state matches the state used in the training
test_set_percent = validation_frac
val_df = normalized_compendium.sample(frac=test_set_percent, random_state=random_state)
val_samples = list(val_df.index)


# In[6]:


# UMAP embedding of original input data

# Get and save model
model = umap.UMAP(random_state=random_state).fit(normalized_compendium)

input_data_UMAPencoded = model.transform(normalized_compendium)
input_data_UMAPencoded_df = pd.DataFrame(data=input_data_UMAPencoded,
                                         index=normalized_compendium.index,
                                         columns=['1','2'])
# Add label
input_data_UMAPencoded_df['dataset'] = 'training'
input_data_UMAPencoded_df.loc[val_samples,'dataset'] = 'validation'

input_data_UMAPencoded_df


# In[12]:


# Plot
fig = ggplot(input_data_UMAPencoded_df, aes(x='1', y='2'))
fig += geom_point(aes(color='dataset'), alpha=0.2)
fig += labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'UMAP of normalized compendium')
fig += theme_bw()
fig += theme(
    legend_title_align = "center",
    plot_background=element_rect(fill='white'),
    legend_key=element_rect(fill='white', colour='white'), 
    legend_title=element_text(family='sans-serif', size=15),
    legend_text=element_text(family='sans-serif', size=12),
    plot_title=element_text(family='sans-serif', size=15),
    axis_text=element_text(family='sans-serif', size=12),
    axis_title=element_text(family='sans-serif', size=15)
    )
fig += guides(colour=guide_legend(override_aes={'alpha': 1}))
fig += scale_color_manual(['#ff6666', '#add8e6'])

print(fig)


# **Observations:**
# * There looks to be a good amount of variance in the compendium overall.
# * Using a split of 25% seems to get a similar distribution of data between training and validation sets.
# * Remember, the dataset is in 17K dimensional space, which will make the small clusters difficult to represent during training
# 
# Overall, having so many features in our dataset, points to the need for more samples to represent the structure in the compendium. For now, we are limited by memory to only select a subset of recount2, but in a future iteration perhaps this will be updated.
