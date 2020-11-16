
# coding: utf-8

# # Coverage of MultiPLIER LV
# 
# The goal of this notebook is to examine why genes were found to be generic. Specifically, this notebook is trying to answer the question: Are generic genes found in more multiplier latent variables compared to specific genes?
# 
# The PLIER model performs a matrix factorization of gene expression data to get two matrices: loadings (Z) and latent matrix (B). The loadings (Z) are constrained to aligned with curated pathways and gene sets specified by prior knowledge [Figure 1B of Taroni et. al.](https://www.cell.com/cell-systems/pdfExtended/S2405-4712\(19\)30119-X). This ensure that some but not all latent variables capture known biology. The way PLIER does this is by applying a penalty such that the individual latent variables represent a few gene sets in order to make the latent variables more interpretable. Ideally there would be one latent variable associated with one gene set unambiguously.
# 
# While the PLIER model was trained on specific datasets, MultiPLIER extended this approach to all of recount2, where the latent variables should correspond to specific pathways or gene sets of interest. Therefore, we will look at the coverage of generic genes versus other genes across these MultiPLIER latent variables, which represent biological patterns.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import random
import textwrap
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ponyo import utils
from generic_expression_patterns_modules import process


# In[2]:


# Get data directory containing gene summary data
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
data_dir = os.path.join(base_dir, "human_general_analysis")

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]


# In[3]:


# Output file
nonzero_figure_filename = "nonzero_LV_coverage.svg"
highweight_figure_filename = "highweight_LV_coverage.svg"


# ## Load data

# In[4]:


# Get gene summary file
ls_data_filename = process.get_gene_summary_file(data_dir)

# Check that there is only one file returned
assert len(ls_data_filename) == 1


# In[5]:


# Load gene summary data
data = pd.read_csv(ls_data_filename[0], sep="\t", index_col=0, header=0)

# Check that genes are unique since we will be using them as dictionary keys below
assert(data.shape[0] == len(data["Gene ID"].unique()))


# In[6]:


# Load multiplier models
# Converted formatted pickle files (loaded using phenoplier environment) from
# https://github.com/greenelab/phenoplier/blob/master/nbs/01_preprocessing/005-multiplier_recount2_models.ipynb
# into .tsv files
multiplier_model_z = pd.read_csv("multiplier_model_z.tsv", sep="\t", index_col=0, header=0)


# In[7]:


# Get a rough sense for how many genes contribute to a given LV
# (i.e. how many genes have a value > 0 for a given LV)
(multiplier_model_z > 0).sum().sort_values(ascending=True)


# ## Get gene data
# 
# Define generic genes based on simulated gene ranking. Refer to [figure](https://github.com/greenelab/generic-expression-patterns/blob/master/human_general_analysis/gene_ranking_log2FoldChange.svg) as a guide.
# 
# **Definitions:**
# * Generic genes: `Rank (simulated) >= 10000` 
# 
# (Having a high rank indicates that these genes are consistently changed across simulated experiments.)
# 
# * Other genes: `Rank (simulated) < 10000` 
# 
# (Having a lower rank indicates that these genes are not consistently changed across simulated experiments - i.e. the genes are specifically changed in an experiment. It could also indicate genes that are consistently unchanged.)

# In[8]:


generic_threshold = 10000
dict_genes = process.get_generic_specific_genes(data, generic_threshold)


# In[9]:


# Check overlap between multiplier genes and our genes
multiplier_genes = list(multiplier_model_z.index)
our_genes = list(data.index)
shared_genes = set(our_genes).intersection(multiplier_genes)

print(len(our_genes))
print(len(shared_genes))


# In[10]:


# Drop gene ids not used in multiplier analysis
processed_dict_genes = process.process_generic_specific_gene_lists(dict_genes, multiplier_model_z)


# In[11]:


# Check numbers add up
assert len(shared_genes) == len(processed_dict_genes["generic"]) + len(processed_dict_genes["other"])


# ## Get coverage of LVs
# 
# For each gene (generic or other) we want to find:
# 1. The number of LVs that gene is present
# 2. The number of LVs that the gene contributes a lot to (i.e. the gene is highly weighted within that LV)

# ### Nonzero LV coverage

# In[12]:


dict_nonzero_coverage = process.get_nonzero_LV_coverage(processed_dict_genes, multiplier_model_z)


# In[13]:


# Check genes mapped correctly
assert processed_dict_genes["generic"][0] in dict_nonzero_coverage["generic"].index
assert len(dict_nonzero_coverage["generic"]) == len(processed_dict_genes["generic"])
assert len(dict_nonzero_coverage["other"]) == len(processed_dict_genes["other"])


# ### High weight LV coverage

# In[14]:


# Quick look at the distribution of gene weights per LV
sns.distplot(multiplier_model_z["LV2"])
plt.yscale("log")


# In[15]:


dict_highweight_coverage = process.get_highweight_LV_coverage(processed_dict_genes, multiplier_model_z, 0.9)


# In[16]:


# Check genes mapped correctly
assert processed_dict_genes["generic"][0] in dict_highweight_coverage["generic"].index
assert len(dict_highweight_coverage["generic"]) == len(processed_dict_genes["generic"])
assert len(dict_highweight_coverage["other"]) == len(processed_dict_genes["other"])


# ### Assemble LV coverage and plot

# In[17]:


all_coverage = []
for gene_label in dict_genes.keys():
    merged_df = pd.DataFrame(
        dict_nonzero_coverage[gene_label],
        columns= ["nonzero LV coverage"]
    ).merge(
        pd.DataFrame(
            dict_highweight_coverage[gene_label],
            columns= ["highweight LV coverage"]
        ),
        left_index=True,
        right_index=True
    ) 
    merged_df['gene type'] = gene_label
    all_coverage.append(merged_df)

all_coverage_df = pd.concat(all_coverage)


# In[18]:


all_coverage_df.head()


# In[19]:


# Plot coverage distribution given list of generic coverage, specific coverage
nonzero_fig = sns.boxplot(
    data=all_coverage_df, 
    x='gene type',
    y='nonzero LV coverage',
    notch=True,
    palette=['powderblue', 'grey']
                         )
plt.ylim(0, 700)
nonzero_fig.set_xlabel("Gene Type",fontsize=14)
nonzero_fig.set_ylabel(textwrap.fill("Number of LVs", width=30),fontsize=14)
nonzero_fig.tick_params(labelsize=14)
nonzero_fig.set_title("Number of LVs genes are present in", fontsize=16)


# In[20]:


# Plot coverage distribution given list of generic coverage, specific coverage
highweight_fig = sns.boxplot(data=all_coverage_df, 
                             x='gene type',
                             y='highweight LV coverage',
                             notch=True,
                             palette=['powderblue', 'grey']
                            )
plt.ylim(0, 700)
highweight_fig.set_xlabel("Gene Type",fontsize=14)
highweight_fig.set_ylabel(textwrap.fill("Number of LVs", width=30),fontsize=14)
highweight_fig.tick_params(labelsize=14)
highweight_fig.set_title("Number of LVs genes contribute highly to", fontsize=16)


# ## Get LVs that generic genes are highly weighted in
# 
# Currently, I get high weight genes if the gene weight > quantile threshold for that specific LV. Each LV has a different threshold based on its own distribution. So each LV has the same number of high weight genes. 
# 
# Should I have normalized the data per LV before calculating which genes are high weight per LV?
# 
# If I normalized then I could probably find a subset of LV that have high weight generic genes. Without this normalization all LVs have some high weight generic genes.

# In[21]:


#thresholds_per_LV = multiplier_model_z.quantile(0.9)
#gene_ids = processed_dict_genes["generic"]
#multiplier_model_z[(multiplier_model_z > thresholds_per_LV)["LV1"] == True]


# In[22]:


#thresholds_per_LV


# ## Save

# In[23]:


# Save plot
nonzero_fig.figure.savefig(
        nonzero_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )

# Save plot
highweight_fig.figure.savefig(
        highweight_figure_filename,
        format="svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0,
        dpi=300,
    )


# **Takeaway:**
# * Generic and other genes have are present in a similar number of LVs. This isn't surprising since the number of genes that contribute to each LV is <1000.
# * Other genes are highly weighted in more LVs compared to generic genes
# * So, generic genes contribute a little to many LVs versus other genes that contribute a lot to some LVs
# * The LVs that were found to contribute alot to can be found in [table](generic_only_LV_summary.tsv). These LVs include -------mainly immune response pathways (monocytes, mast cell activation), wound healing (collagen formation), cell signaling (focal adhesion, integrin1) 
# 
# **Overall, it looks like generic genes are associated with many pathways, acting as *gene hubs*, which is why they are "generic"**
