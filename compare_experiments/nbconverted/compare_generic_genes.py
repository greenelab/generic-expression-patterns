
# coding: utf-8

# # Compare generic genes
# 
# The goal of this notebook is to compare the generic genes found using 2 different recount2 template experiments to determine ____

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import seaborn as sns
import pandas as pd
from ponyo import utils


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]

project_id1 = "SRP012656"
project_id2 = "SRP061689"


# In[3]:


# Get data directory containing gene summary data
data_dir = os.path.join(base_dir, "human_general_analysis")

# Get gene ranking files
gene_ranking_filename1 = os.path.join(data_dir, f"generic_gene_summary_{project_id1}.tsv")
gene_ranking_filename2 = os.path.join(data_dir, f"generic_gene_summary_{project_id2}.tsv")

# Get template data
template_filename1 = os.path.join(data_dir, "data", f"processed_recount2_template_{project_id1}.tsv")
template_filename2 = os.path.join(data_dir, "data", f"processed_recount2_template_{project_id2}.tsv")


# ## Correlation between rankings

# In[4]:


# Load gene ranking
gene_ranking_summary1 = pd.read_csv(gene_ranking_filename1, sep="\t", index_col=0, header=0)
gene_ranking_summary2 = pd.read_csv(gene_ranking_filename2, sep="\t", index_col=0, header=0)


# In[5]:


# Get simulated ranking
gene_ranking1 = gene_ranking_summary1["Rank (simulated)"].rename("Rank 1")
gene_ranking2 = gene_ranking_summary2["Rank (simulated)"].rename("Rank 2")

# Combine ranking
gene_ranking_combined = pd.concat([gene_ranking1, gene_ranking2], axis=1)


# In[16]:


gene_ranking_combined.loc[gene_ranking_combined.isna().any(axis=1)]


# In[7]:


# Plot correlation between ranking
fig = sns.jointplot(
    data=gene_ranking_combined,
    x="Rank 1",
    y="Rank 2",
    kind="hex",
    marginal_kws={"color": "white"},
)

fig.set_axis_labels(
    f"Ranking in {project_id1}", f"Ranking in {project_id2}", fontsize=14, fontname="Verdana"
)

"""fig.savefig(
    output_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)"""


# In[8]:


print(gene_ranking1.shape)
gene_ranking_summary1.loc["A1CF"]


# In[9]:


print(gene_ranking2.shape)
gene_ranking_summary2.loc["A1CF"]


# **Takeaway:**
# 
# * Looks like there is good concordance between highly ranked genes (i.e. generic genes)
# * Are the gene rankings robust? We expect that there are some genes that are generally generic and there are also genes that are generic within the given context.

# ## Examine gene expression data

# In[10]:


template_1 = pd.read_csv(template_filename1, sep="\t", index_col=0, header=0)
template_2 = pd.read_csv(template_filename2, sep="\t", index_col=0, header=0)


# In[17]:


sns.distplot(template_1.mean(), kde=False)
sns.distplot(template_2.mean(), kde=False)


# **Takeaway:**
# 
# Are there differences in the distribution of mean gene expression? Doesn't look like there are
# 
# Do we expect that differences in gene expression in the template experiment will yeild different generic genes? Theoretically, if a gene is lowly expressed in the context of template experiment 1 it shouldn't be found to be generic. But this same gene could be found to be generic in the context of template experiment 2 if it is more expressed.
