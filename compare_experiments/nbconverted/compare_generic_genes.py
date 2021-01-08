
# coding: utf-8

# # Compare generic genes
# 
# The goal of this notebook is to compare the generic genes found using 2 different recount2 template experiments.

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


# In[6]:


print(gene_ranking_combined.shape)
gene_ranking_combined.head()


# In[7]:


# Check for NAs
gene_ranking_combined[pd.isnull(gene_ranking_combined).any(axis=1)]


# In[8]:


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

output_figure_filename = "concordance_between_recount2_templates.svg"
fig.savefig(
    output_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)


# **Takeaway:**
# 
# * Looks like there is good concordance between highly ranked genes (i.e. generic genes)
# * In general, we expect that there are some genes that are generally generic (i.e. those that are concordant between these two experiments) and there are also genes that are generic within the given context of the specific experiment.

# ## Examine gene expression data

# In[9]:


# Read expression data
template_1 = pd.read_csv(template_filename1, sep="\t", index_col=0, header=0)
template_2 = pd.read_csv(template_filename2, sep="\t", index_col=0, header=0)


# In[10]:


# Get concordance genes
concordant_genes = list(gene_ranking_combined[(gene_ranking_combined["Rank 1"]>15000) &
                                         (gene_ranking_combined["Rank 2"]>15000)].index)

# Get disconcordant genes
discordant_genes = set(gene_ranking_combined.index).difference(concordant_genes)


# In[11]:


# Distribution of concordant genes in template experiment 1
template1_mean = template_1.mean()

print("Percent concordant genes with 0 expression in template 1:", 
      len(template1_mean[concordant_genes].loc[template1_mean[concordant_genes]==0])/len(template1_mean[concordant_genes]))

print("Percent nonzero concordant genes in template 1:", 
      len(template1_mean[concordant_genes].loc[(template1_mean[concordant_genes]>0) &
                                           (template1_mean[concordant_genes]<1000)])/len(template1_mean[concordant_genes]))
sns.distplot(template_1.mean()[concordant_genes], kde=False)


# In[12]:


# Distribution of concordant genes in template experiment 2
template2_mean = template_2.mean()
print("Percent concordant genes with 0 expression in template 2:",
      len(template2_mean[concordant_genes].loc[template2_mean[concordant_genes]==0])/len(template2_mean[concordant_genes]))

print("Percent nonzero concordant genes in template 2:",
      len(template2_mean[concordant_genes].loc[(template2_mean[concordant_genes]>0) &
                                                (template2_mean[concordant_genes]<1000)])/len(template2_mean[concordant_genes]))

# There are more 0 expressed genes in this template experiment
sns.distplot(template_2.mean()[concordant_genes], kde=False)


# In[13]:


# Distribution of discordant gense in template experiment 1
template1_mean = template_1.mean()

print("Percent discordant genes with 0 expression in template 1:", 
      len(template1_mean[discordant_genes].loc[template1_mean[discordant_genes]==0])/len(template1_mean[discordant_genes]))

print("Percent nonzero discordant genes in template 1:", 
      len(template1_mean[discordant_genes].loc[(template1_mean[discordant_genes]>0) &
                                           (template1_mean[discordant_genes]<1000)])/len(template1_mean[discordant_genes]))

print(len(template1_mean[discordant_genes].loc[template1_mean[discordant_genes]>0])/len(template1_mean[discordant_genes]))
sns.distplot(template_1.mean()[discordant_genes], kde=False)


# In[14]:


# Distribution of discordant genes in template experiment 2
template2_mean = template_2.mean()

print("Percent discordant genes with 0 expression in template 2:",
      len(template2_mean[discordant_genes].loc[template2_mean[discordant_genes]==0])/len(template2_mean[discordant_genes]))

print("Percent nonzero discordant genes in template 2:",
      len(template2_mean[discordant_genes].loc[(template2_mean[discordant_genes]>0) &
                                                (template2_mean[discordant_genes]<1000)])/len(template2_mean[discordant_genes]))

sns.distplot(template_2.mean()[discordant_genes], kde=False)


# **Takeaway:**
# 
# Doesn't appear to be much of a difference between the distribution of average gene expression values for these two experiments. 
# 
# Theoretically, I would expect the scenario where a gene is lowly expressed in the context of template experiment 1 and therefore not found to be generic. But this same gene could be found to be generic in the context of template experiment 2 if it is more expressed. Its possible that differences in gene expression distribution can change which genes are found to be generic given that the simulation is producing experiments with a similar context. 
# 
# In this case, despite having similar gene expression distributions there are still many differences in gene ranking. This suggests to me that level of gene expression activity doesn't matter as much as the overall patterns perhaps. 
