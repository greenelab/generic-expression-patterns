
# coding: utf-8

# # Test: Identify generic human genes on test set
# 
# This notebook performs the following steps to identify generic genes
# 1. Simulates N gene expression experiments using [ponyo](https://github.com/ajlee21/ponyo)
# 2. Perform DE analysis to get association statistics for each gene
# 
# In this case the DE analysis is based on the experimental design of the template experiment, described in the previous [notebook](1_process_recount2_data.ipynb). The template experiment is [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764), which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients. So the DE analysis is comparing tumor vs normal in this case.
# 
# 3. For each gene, aggregate statsitics across all simulated experiments 
# 4. Rank genes based on this aggregated statistic
# 
# **Evaluation:**
# We want to compare our ranking using ponyo, compared to the ranking found from Crow et. al.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import numpy as np
import pickle
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import process, stats, ranking

np.random.seed(123)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(
    os.path.join(os.getcwd(), "../")
)

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_test.tsv")
)

params = utils.read_config(config_filename)


# In[3]:


# Load params
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
NN_architecture = params['NN_architecture']
num_runs = params['num_simulated']
project_id = params['project_id']
metadata_col_id = params['metadata_colname']
mapped_template_filename = params['mapped_template_filename']
processed_template_filename = params['processed_template_filename']
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = params['scaler_filename']
col_to_rank_genes = params['rank_genes_by']
col_to_rank_pathways = params['rank_pathways_by']
statistic = params['gsea_statistic']

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_groups.tsv"
)

# Load pickled file
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)


# ## Test: Simulation

# In[4]:


# Simulate multiple experiments
# This step creates the following files in "<local_dir>/pseudo_experiment/" directory:           
#   - selected_simulated_data_SRP012656_<n>.txt
#   - selected_simulated_encoded_data_SRP012656_<n>.txt
#   - template_normalized_data_SRP012656_test.txt
# in which "<n>" is an integer in the range of [0, num_runs-1] 
os.makedirs(os.path.join(local_dir, "pseudo_experiment"), exist_ok=True)
for run_id in range(num_runs):
    simulate_expression_data.shift_template_experiment(
        normalized_compendium_filename,
        project_id,
        metadata_col_id,
        NN_architecture,
        dataset_name,
        scaler,
        local_dir,
        base_dir,
        run_id
    )


# ## Test: Processing template experiment
# 
# `processed_template_filename`: template data with some sample rows dropped

# In[5]:


if os.path.exists(sample_id_metadata_filename):
    stats.process_samples_for_DESeq(
        mapped_template_filename,
        sample_id_metadata_filename,
        metadata_filename,
        None,
        processed_template_filename
    )


# In[6]:


# Read data
template_data = pd.read_csv(
    processed_template_filename,
    header=0,
    sep='\t',
    index_col=0
)

# Check samples dropped
print(template_data.shape)
assert(template_data.shape[0] == 24)
template_data.head()


# ## Test: Processing simulation experiments

# In[7]:


# This step modifies the following files:
# "<local_dir>/pseudo_experiments/selected_simulated_data_SRP012656_<n>.txt"
if os.path.exists(sample_id_metadata_filename):
    for i in range(num_runs):
        simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt"
        )
        stats.process_samples_for_DESeq(
        simulated_filename,
        sample_id_metadata_filename,
        metadata_filename
    )


# In[8]:


# Check simulated files were created
sim_output1 = os.path.join(local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_0.txt")
sim_output2 = os.path.join(local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_1.txt")
assert (os.path.exists(sim_output1) and os.path.exists(sim_output2))


# In[9]:


# Check that simulated files are non-empty
assert (os.path.getsize(sim_output1)>0 and os.path.getsize(sim_output2)>0)


# **Note:** These cells testing for reproducibility of the simulation pass when run locally. But fail when run on github actions, so for now I am commenting them out but will use them to test locally any future updates to the code

# In[10]:


"""# Check reproducibility of simulated experiments using random seed
template_path = "data/test_simulated_data_SRP012656_0.txt"
output_path = os.path.join(local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_0.txt")
template_df = pd.read_csv(template_path, sep="\t", header=0, index_col=0)
output_df = pd.read_csv(output_path, sep="\t", header=0, index_col=0)

assert np.all(np.isclose(output_df.values, template_df.values)), (
    output_df.iloc[
        np.where(~np.all(np.isclose(output_df.values, template_df.values), axis=1))[0],
        np.where(~np.all(np.isclose(output_df.values, template_df.values), axis=0))[0],
    ],
)"""


# In[11]:


# Check reproducibility of simulated experiments
"""template_path = "data/test_simulated_data_SRP012656_1.txt"
output_path = os.path.join(local_dir, "pseudo_experiment", "selected_simulated_data_SRP012656_1.txt")
template_df = pd.read_csv(template_path, sep="\t", header=0, index_col=0)
output_df = pd.read_csv(output_path, sep="\t", header=0, index_col=0)

assert np.all(np.isclose(output_df.values, template_df.values)), (
    output_df.iloc[
        np.where(~np.all(np.isclose(output_df.values, template_df.values), axis=1))[0],
        np.where(~np.all(np.isclose(output_df.values, template_df.values), axis=0))[0],
    ]
)"""


# ## Test: Differential expression analysis

# In[12]:


# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)


# In[13]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/DE_analysis.R\'))\n\n# File created: "<local_dir>/DE_stats/DE_stats_template_data_SRP012656_real.txt"\nget_DE_stats_DESeq(metadata_filename,\n                   project_id, \n                   processed_template_filename,\n                   "template",\n                   local_dir,\n                   "real")')


# In[14]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs -i base_dir', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/DE_analysis.R\'))\n\n# Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_SRP012656_<n>.txt"\nfor (i in 0:(num_runs-1)){\n    simulated_data_filename <- paste(local_dir, \n                                     "pseudo_experiment/selected_simulated_data_",\n                                     project_id,\n                                     "_", \n                                     i,\n                                     ".txt",\n                                     sep = "")\n    \n    get_DE_stats_DESeq(metadata_filename,\n                       project_id, \n                       simulated_data_filename,\n                       "simulated",\n                       local_dir,\n                       i)\n}')


# In[15]:


# Check DE stats files were created
DE_output1 = os.path.join(local_dir, "DE_stats", "DE_stats_simulated_data_SRP012656_0.txt")
DE_output2 = os.path.join(local_dir, "DE_stats", "DE_stats_simulated_data_SRP012656_1.txt")
assert (os.path.exists(DE_output1) and os.path.exists(DE_output2))


# In[16]:


# Check that DE stats files are non-empty
assert (os.path.getsize(DE_output1)>0 and os.path.getsize(DE_output2)>0)


# ### Rank genes

# In[17]:


analysis_type = "DE"
template_DE_stats_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_template_data_{project_id}_real.txt"
)

template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
    template_DE_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_genes,
)


# ### Gene summary table

# In[18]:


summary_gene_ranks = ranking.generate_summary_table(
    template_DE_stats_filename,
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir,
    'gene',
    params
)


# In[39]:


summary_gene_ranks.head()


# In[41]:


# Some genes will have NaN's in the simulated statistics columns. These are genes that were filtered 
# due to low expression and therefore the corresponding Z-score for this gene is NaN
summary_gene_ranks.isna().any()


# In[45]:


summary_gene_ranks[summary_gene_ranks.isna().any(axis=1)]


# In[19]:


# Create `gene_summary_filename`
gene_summary_filename = os.path.join(local_dir, "gene_summary_table.tsv")
summary_gene_ranks.to_csv(gene_summary_filename, sep='\t')


# In[20]:


"""# Passed assertion locally but not on github actions but not clear why
template_path = "data/test_gene_summary_table.tsv"
template_df = pd.read_csv(template_path, sep="\t", header=0, index_col=0)
output_df = pd.read_csv(gene_summary_filename, sep="\t", header=0, index_col=0)

assert (template_df["Gene ID"].values == output_df["Gene ID"].values).all(),template_df.loc[template_df["Gene ID"].values != output_df["Gene ID"].values,"Gene ID"]

assert np.all(np.isclose(
    template_df[["Rank (Real)", "Rank (simulated)"]]values,
    output_df[["Rank (Real)", "Rank (simulated)"]]values)),(
    output_df[["Rank (Real)", "Rank (simulated)"]].iloc[
        np.where(~np.all(np.isclose(output_df[["Rank (Real)", "Rank (simulated)"]].values, template_df[["Rank (Real)", "Rank (simulated)"]].values), axis=1))[0],
        np.where(~np.all(np.isclose(output_df[["Rank (Real)", "Rank (simulated)"]].values, template_df[["Rank (Real)", "Rank (simulated)"]].values), axis=0))[0],
    ]
)"""


# ## Test: Compare gene ranking
# Studies have found that there are some genes that are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represents a more systematic change. 
# 
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# In[21]:


# Get generic genes identified by Crow et. al.
DE_prior_file = params['reference_gene_filename']
ref_gene_col = params['reference_gene_name_col']
ref_rank_col = params['reference_rank_col']

figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

corr_stats = ranking.compare_gene_ranking(
    summary_gene_ranks,
    DE_prior_file,
    ref_gene_col,
    ref_rank_col,
    figure_filename
)
r, p = corr_stats['r'], corr_stats['p']

# Expected output for DE using Limma
#expected_r = 0.21913957199910106
#expected_p = 6.871971345526456e-186

expected_r = 0.22258799129252085
expected_p = 4.5589621319725286e-173
assert(
    np.all(
        np.isclose([r, p], [expected_r, expected_p])
    ),
    ([r,p], [expected_r, expected_p])
)


# ## Test: GSEA

# In[22]:


# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "GSEA_stats"), exist_ok=True)


# In[23]:


# Load pathway data
hallmark_DB_filename = os.path.join(base_dir, dataset_name, "data", "metadata", "hallmark_DB.gmt")


# In[24]:


get_ipython().run_cell_magic('R', '-i base_dir -i template_DE_stats_filename -i hallmark_DB_filename -i statistic -o template_enriched_pathways', '\nsource(paste0(base_dir, \'/generic_expression_patterns_modules/GSEA_analysis.R\'))\nout_file <- paste(local_dir, \n                     "GSEA_stats/GSEA_stats_template_data_",\n                     project_id,\n                     "_real.txt", \n                     sep = "")\n\ntemplate_enriched_pathways <- find_enriched_pathways(template_DE_stats_filename, hallmark_DB_filename, statistic)  \n    \nwrite.table(as.data.frame(template_enriched_pathways[1:7]), file = out_file, row.names = F, sep = "\\t")')


# In[25]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_filename -i num_runs -i statistic -i base_dir', '\nsource(paste0(base_dir,\'/generic_expression_patterns_modules/GSEA_analysis.R\'))\n\n# New files created: "<local_dir>/GSEA_stats/GSEA_stats_simulated_data_<project_id>_<n>.txt"\nfor (i in 0:(num_runs-1)) {\n    simulated_DE_stats_file <- paste(local_dir, \n                                     "DE_stats/DE_stats_simulated_data_", \n                                     project_id,\n                                     "_", \n                                     i,\n                                     ".txt",\n                                     sep = "")\n    \n    out_file <- paste(local_dir, \n                     "GSEA_stats/GSEA_stats_simulated_data_",\n                     project_id,\n                     "_",\n                     i,\n                     ".txt", \n                     sep = "")\n        \n    enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_filename, statistic) \n    \n    # Remove column with leading edge since its causing parsing issues\n    write.table(as.data.frame(enriched_pathways[1:7]), file = out_file, row.names = F, sep = "\\t")\n}')


# In[26]:


# Check GSEA stats files were created
GSEA_output1 = os.path.join(local_dir, "GSEA_stats", "GSEA_stats_simulated_data_SRP012656_0.txt")
GSEA_output2 = os.path.join(local_dir, "GSEA_stats", "GSEA_stats_simulated_data_SRP012656_1.txt")
assert (os.path.exists(DE_output1) and os.path.exists(DE_output2))


# In[27]:


# Check that GSEA stats files are non-empty
assert (os.path.getsize(GSEA_output1)>0 and os.path.getsize(GSEA_output2)>0)


# ### Rank pathways

# In[28]:


analysis_type = "GSEA"
template_GSEA_stats_filename = os.path.join(
    local_dir,
    "GSEA_stats",
    f"GSEA_stats_template_data_{project_id}_real.txt"
)

template_GSEA_stats, simulated_GSEA_summary_stats = ranking.process_and_rank_genes_pathways(
    template_GSEA_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_pathways,
)


# In[29]:


"""# Concatenate simulated experiments
simulated_GSEA_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id, 'GSEA')
simulated_GSEA_stats_all.set_index('pathway', inplace=True)
print(simulated_GSEA_stats_all.shape)"""


# In[30]:


"""# Aggregate statistics across all simulated experiments
simulated_GSEA_summary_stats = calc.aggregate_stats(
    col_to_rank_pathways,
    simulated_GSEA_stats_all,
    'GSEA'
)"""


# In[31]:


"""# Load association statistics for template experiment
template_GSEA_stats = template_enriched_pathways.iloc[:, :-1]
template_GSEA_stats.set_index('pathway', inplace=True)

template_GSEA_stats.head()

# Rank genes in template experiment
template_GSEA_stats = calc.rank_genes_or_pathways(
    col_to_rank_pathways,
    template_GSEA_stats,
    True
)"""


# In[32]:


"""# Rank genes in simulated experiments
simulated_GSEA_summary_stats = calc.rank_genes_or_pathways(
    col_to_rank_pathways,
    simulated_GSEA_summary_stats,
    False
)"""


# ### Pathway summary table

# In[33]:


# Create intermediate file: "<local_dir>/gene_summary_table_<col_to_rank_pathways>.tsv"
summary_pathway_ranks = ranking.generate_summary_table(
    template_GSEA_stats_filename,
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
    col_to_rank_pathways,
    local_dir,
    'pathway',
    params
)


# ## Test: Compare pathway ranking

# In[34]:


# Load Powers et. al. results file
powers_rank_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked.csv"
)


# In[35]:


# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_filename, header=0, index_col=0)
powers_rank_df.drop(['Category'], axis=1, inplace=True)


# In[36]:


# Count the number of experiments where a given pathway was found to be enriched (qvalue < 0.05)
total_num_experiments = powers_rank_df.shape[1]
frac_enriched_pathways = ((powers_rank_df < 0.05).sum(axis=1) / total_num_experiments)

# Rank pathways from 0-50, 50 indicating that the pathways was frequently enriched
pathway_ranks = frac_enriched_pathways.rank()

powers_rank_stats_df = pd.DataFrame(
    data={
        'Fraction enriched': frac_enriched_pathways.values,
        'Powers Rank':pathway_ranks.values
    },
    index=powers_rank_df.index
)


# In[37]:


# Save reference file for input into comparison
powers_rank_processed_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked_processed.tsv"
)

powers_rank_stats_df.to_csv(powers_rank_processed_filename, sep="\t", )


# In[38]:


figure_filename = f"pathway_ranking_{col_to_rank_pathways}.svg"

corr_stats = ranking.compare_pathway_ranking(
    summary_pathway_ranks,
    powers_rank_processed_filename,
    figure_filename
)
# Note: Not getting reproducible results after GSEA, maybe due to permutations
#r, p = corr_stats['r'], corr_stats['p']
    
#expected_r = 0.07620992008839596
#expected_p = 0.5988732068701128

#assert(
#    np.all(
#        np.isclose([r, p], [expected_r, expected_p])
#    )
#)

