
# coding: utf-8

# # Process MULTIPLIER model data
# 
# Raw multiplier model data can be found [here](https://github.com/greenelab/multi-plier).
# 
# This multiplier model (Robject) was formatted in python in [phenoplier repo](https://github.com/greenelab/phenoplier/blob/master/nbs/01_preprocessing/005-multiplier_recount2_models.ipynb). The python pickled data files were output to a [shared google drive](https://drive.google.com/drive/folders/12wvqGzFpFBUOX_CsFkFvvcbKfhl648LE). Since these pickle files were generated using a different conda environment, this notebook is loading the data from the pickle files into .tsv files for use in our analysis.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import pickle

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


# In[2]:


multiplier_dir = os.path.join("/home/alexandra/Documents/Data/Generic_expression_patterns/", "multiplier")


# # Load pickle files and save as .tsv

# In[3]:


multiplier_u_filename = os.path.join(multiplier_dir, "multiplier_model_u.pkl")
multiplier_z_filename = os.path.join(multiplier_dir, "multiplier_model_z.pkl")


# In[4]:


# Load pickled file
with open(multiplier_u_filename, "rb") as scaler_fh:
    multiplier_model_u = pickle.load(scaler_fh)

with open(multiplier_z_filename, "rb") as scaler_fh:
    multiplier_model_z = pickle.load(scaler_fh)


# In[5]:


print(multiplier_model_u.shape)
multiplier_model_u.head()


# In[6]:


print(multiplier_model_z.shape)
multiplier_model_z.head()


# In[7]:


# make sure I'm seeing the same when loaded with R
assert multiplier_model_z.loc['GAS6', 'LV2'] == 0
assert multiplier_model_z.loc['GAS6', 'LV3'] == 0.039437739697954444
assert multiplier_model_z.loc['GAS6', 'LV984'] == 0.3473620915326928
assert multiplier_model_z.loc['GAS6', 'LV987'] == 0

assert multiplier_model_z.loc['SPARC', 'LV981'] == 0
assert multiplier_model_z.loc['SPARC', 'LV986'].round(8) == 0.12241734


# In[8]:


# Save
multiplier_model_u.to_csv("multiplier_model_u.tsv", sep="\t")
multiplier_model_z.to_csv("multiplier_model_z.tsv", sep="\t")


# # Format multiplier summary data

# In[9]:


readRDS = ro.r['readRDS']


# In[10]:


multiplier_full_model = readRDS(os.path.join(multiplier_dir,
                                             "recount_PLIER_model.RDS"))


# In[11]:


multiplier_model_matrix = multiplier_full_model.rx2('summary')


# In[12]:


with localconverter(ro.default_converter + pandas2ri.converter):
  multiplier_model_matrix_values = ro.conversion.rpy2py(multiplier_model_matrix)


# In[13]:


multiplier_model_matrix_df = pd.DataFrame(
    data=multiplier_model_matrix_values,
    index=multiplier_model_matrix.rownames,
    columns=multiplier_model_matrix.colnames
)


# In[14]:


multiplier_model_matrix_df.head()


# In[15]:


# Save
multiplier_model_matrix_df.to_csv("multiplier_model_summary.tsv", sep="\t")

