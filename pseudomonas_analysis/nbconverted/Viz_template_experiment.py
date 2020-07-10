
# coding: utf-8

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
                      geom_point,
                      aes, 
                      ggsave, 
                      theme_bw,
                      theme,
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

np.random.seed(123)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))
config_file = os.path.abspath(os.path.join(base_dir,
                                           "config_pseudomonas.tsv"))
params = utils.read_config(config_file)


# In[3]:


# Load parameters
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
project_id = params['project_id']
template_data_file = params['template_data_file']


# In[4]:


# Load metadata file with grouping assignments for samples
metadata_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    project_id+"_groups.tsv")


# In[5]:


# Read template data
data = pd.read_csv(template_data_file, sep="\t", header=0, index_col=0)

data.head()


# In[6]:


# Read metadata
metadata = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)

metadata.head()


# In[7]:


# Embed expression data into low dimensional space
model = umap.UMAP(random_state=123).fit(data)
data_encoded = model.transform(data)

data_encoded_df = pd.DataFrame(data=data_encoded,
                               index=data.index,
                               columns=['1','2'])


# In[8]:


# Label samples
group1_ids = list(metadata[metadata['group']==1].index)

#data_encoded_df['group'] = 'clinical multi-drug resistant'
#data_encoded_df.loc[group1_ids,'group'] = 'clinical'
#data_encoded_df.loc['GSM625982.CEL','group'] = 'control'
data_encoded_df['group'] = 'untreated'
data_encoded_df.loc[group1_ids,'group'] = 'treated with tobramycin'


# In[9]:


data_encoded_df.head()


# In[10]:


# Plot PAO1
fig = ggplot(data_encoded_df, aes(x='1', y='2'))
fig += geom_point(aes(color='group'), alpha=0.7)
fig += labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'Gene expression of template experiment')
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

print(fig)

