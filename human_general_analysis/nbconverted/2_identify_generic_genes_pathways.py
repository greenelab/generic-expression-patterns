#!/usr/bin/env python
# coding: utf-8

# # Identify generic genes and pathways
# 
# **This notebook performs the following steps to identify generic genes:**
# 1. Simulates N gene expression experiments using [ponyo](https://github.com/ajlee21/ponyo)
# 2. Perform DE analysis to get association statistics for each gene
# 
#   In this case the DE analysis is based on the experimental design of the template experiment, described in the previous [notebook](1_process_recount2_data.ipynb). The template experiment is [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764), which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients. So the DE analysis is comparing tumor vs normal in this case.
# 
# 
# 3. For each gene, aggregate statsitics across all simulated experiments 
# 4. Rank genes based on this aggregated statistic (i.e. log fold change, or p-value)
# 
# 
# **This notebook performs the following steps to identify generic gene sets (pathways):**
# 1. Using the same simulated experiments from above, perform GSEA analysis. This analysis will determine whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list of genes, where genes are ranked by log fold change, indicating a correlation with change in expression.
# 2. For each gene set (pathway), aggregate statistics across all simulated experiments
# 3. Rank gene sets based on this aggregated statistic
# 
# **Evaluation:**
# * We want to compare the ranking of genes identified using the above method with the ranking found from Crow et. al., which identified a set of genes as generic based on their ranking;
# * We want to compare the ranking of pathways identified using the above method with the ranking found from Powers et. al., which identified a set of pathways as generic based on their ranking;
# * This comparison will validate our method being used as a way to automatically identify generic genes and pathways.

# In[ ]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import sys
import pandas as pd
import numpy as np
import pickle

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import calc, process

np.random.seed(123)


# In[ ]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human.tsv")
)

params = utils.read_config(config_filename)


# In[ ]:


# Load params
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
NN_architecture = params['NN_architecture']
num_runs = params['num_simulated']
project_id = params['project_id']
metadata_col_id = params['metadata_colname']
processed_template_filename = params['processed_template_filename']
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = params['scaler_filename']
col_to_rank_genes = params['rank_genes_by']
col_to_rank_pathways = params['rank_pathways_by']
compare_genes = params['compare_genes']
statistic = params['gsea_statistic']

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv"
)

# Load pickled file
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)


# In[ ]:


# Output files
gene_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_gene_summary_{project_id}.tsv"
)

pathway_summary_filename = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}.tsv"
)


# ### Simulate experiments using selected template experiment

# In[ ]:


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


# Since this experiment contains both RNA-seq and smRNA-seq samples which are in different ranges, we will drop smRNA samples so that samples are within the same range. The analysis identifying these two subsets of samples can be found in this [notebook](../explore_data/0_explore_input_data.ipynb)

# In[ ]:


# This step modifies the following files:
# "<local_dir>/pseudo_experiments/selected_simulated_data_SRP012656_<n>.txt"
if os.path.exists(sample_id_metadata_filename):
    # Read in metadata
    metadata = pd.read_csv(sample_id_metadata_filename, sep='\t', header=0, index_col=0)
    
    # Get samples to be dropped
    sample_ids_to_drop = list(metadata[metadata["processing"] == "drop"].index)

    process.subset_samples(
        sample_ids_to_drop,
        num_runs,
        local_dir,
        project_id
    )


# In[ ]:


# Round simulated read counts to int in order to run DESeq.
# This step modifies the following files again:
# "<local_dir>/pseudo_experiments/selected_simulated_data_SRP012656_<n>.txt"
process.recast_int(num_runs, local_dir, project_id)


# ### Differential expression analysis

# In[ ]:


# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_groups.tsv"
)


# In[ ]:


# Check whether ordering of sample ids is consistent between gene expression data and metadata
process.compare_and_reorder_samples(processed_template_filename, metadata_filename)


# In[ ]:


# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)


# In[ ]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i processed_template_filename -i local_dir', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\n# File created: "<local_dir>/DE_stats/DE_stats_template_data_SRP012656_real.txt"\nget_DE_stats_DESeq(metadata_filename,\n                   project_id, \n                   processed_template_filename,\n                   "template",\n                   local_dir,\n                   "real")')


# In[ ]:


# Check number of DEGs
template_DE_stats_filename = os.path.join(
    local_dir,
    "DE_stats",
    f"DE_stats_template_data_{project_id}_real.txt"
)

template_DE_stats = pd.read_csv(
    template_DE_stats_filename, 
    sep="\t", 
    header=0, 
    index_col=0
)

selected = template_DE_stats[(template_DE_stats['padj']<0.01) & (abs(template_DE_stats['log2FoldChange'])>1)]
print(selected.shape)


# In[ ]:


# Check whether ordering of sample ids is consistent between gene expression data and metadata
for i in range(num_runs):
    simulated_data_filename = os.path.join(
        local_dir,
        "pseudo_experiment",
        f"selected_simulated_data_{project_id}_{i}.txt"
    )
        
    process.compare_and_reorder_samples(simulated_data_filename, metadata_filename)


# In[ ]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\n# Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_SRP012656_<n>.txt"\nfor (i in 0:(num_runs-1)){\n    simulated_data_filename <- paste(local_dir, \n                                     "pseudo_experiment/selected_simulated_data_",\n                                     project_id,\n                                     "_", \n                                     i,\n                                     ".txt",\n                                     sep = "")\n    \n    get_DE_stats_DESeq(metadata_filename,\n                       project_id, \n                       simulated_data_filename,\n                       "simulated",\n                       local_dir,\n                       i)\n}')


# **Validation:**
# * As a quick validation, [Kim et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3566005/) found 1459 DEGs (543 upregulated and 916 downregulated) using used the Bowtie and NEUMA applications for the mapping and quantification of RNA-Seq data. They used *edgeR* with a rigorous filtering procedure based on false discovery rates, minimum applicable patient numbers, and gene expression levels was devised to select reliable sets of DEGs and DEIs (see File S8 for details). For the
# 
# * Our results found ~3K DEGs which is close enough in range given that the data was processed using different methods. recount2 resource were aligned with the splice-aware Rail-RNA aligner

# ### Rank genes

# In[ ]:


# Concatenate simulated experiments
simulated_DE_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id, 'DE')

print(simulated_DE_stats_all.shape)


# In[ ]:


# Take absolute value of logFC and t statistic
simulated_DE_stats_all = process.abs_value_stats(simulated_DE_stats_all)


# In[ ]:


# Aggregate statistics across all simulated experiments
simulated_DE_summary_stats = calc.aggregate_stats(
    col_to_rank_genes,
    simulated_DE_stats_all,
    'DE'
)


# In[ ]:


# Take absolute value of logFC and t statistic
template_DE_stats = process.abs_value_stats(template_DE_stats)

# Rank genes in template experiment
template_DE_stats = calc.rank_genes_or_pathways(
    col_to_rank_genes,      
    template_DE_stats,
    True
)


# In[ ]:


# Rank genes in simulated experiments
simulated_DE_summary_stats = calc.rank_genes_or_pathways(
    col_to_rank_genes,
    simulated_DE_summary_stats,
    False
)


# ### Gene summary table

# In[ ]:


summary_gene_ranks = process.generate_summary_table(
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir
)

summary_gene_ranks.head()


# In[ ]:


# Create `gene_summary_fielname`
summary_gene_ranks.to_csv(gene_summary_filename, sep='\t')


# ### GSEA 
# **Goal:** To detect modest but coordinated changes in prespecified sets of related genes (i.e. those genes in the same pathway or share the same GO term).
# 
# 1. Rank all genes using DE association statistics. In this case we used the p-value scores to rank genes. logFC returned error -- need to look into this.
# 2. An enrichment score (ES) is defined as the maximum distance from the middle of the ranked list. Thus, the enrichment score indicates whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list (indicating a correlation with change in expression). 
# 3. Estimate the statistical significance of the ES by a phenotypic-based permutation test in order to produce a null distribution for the ES (i.e. scores based on permuted phenotype)

# In[ ]:


# Load pathway data
hallmark_DB_filename = os.path.join(local_dir, "hallmark_DB.gmt")


# In[ ]:


get_ipython().run_cell_magic('R', '-i base_dir -i template_DE_stats_filename -i hallmark_DB_filename -i statistic -o template_enriched_pathways', "\nsource(paste(base_dir, 'generic_expression_patterns_modules/GSEA_analysis.R', sep='/'))\ntemplate_enriched_pathways <- find_enriched_pathways(template_DE_stats_filename, hallmark_DB_filename, statistic)")


# In[ ]:


print(template_enriched_pathways.shape)
template_enriched_pathways[template_enriched_pathways['padj'] < 0.05].sort_values(by='padj').head()


# In[ ]:


template_enriched_pathways[template_enriched_pathways['padj'] >= 0.05].sort_values(by='padj')[:15]


# ### Quick validation
# * GSEA seems to be working as expected when we visualized the distribution of enriched pathways. Distribution of enrichment scores are consistent with what we would expect for not significant pathways (i.e. genes are evenly distributed) and significant pathways (i.e. genes are clustered to one end). 
# 
# <img src="GSEA examination.png" width="800">
# 
# * [Kim et. al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3566005/) *used the Ingenuity Pathway Analysis (IPA) software, which uses a database of evidence manually compiled from the literature. The most enriched term in the diseases and disorders category was cancer (p value = 2.13E-42), which supports the validity of our gene set. Other relevant terms in the molecular and cellular functions category included cellular growth and proliferation (p value = 1.71E-17) and cell death (p value = 1.97E-17). The IPA results are presented in Figure S9 in File S8. Gene ontology (GO) analysis produced similar results to IPA, albeit in a less comprehensive manner (data not shown).*
# 
# * We found the following pathways as being signifcantly enriched in DEGs: HALLMARK_G2M_CHECKPOINT, TNFA_SIGNALING, INFLAMMATORY_RESPONSE, APOPTOSIS, HYPOXIA. These pathways are expected and consistent with what the publication found (cancer pathways, cell death pathways, cell proliferation).
# 
# * However, there are a few pathways that were not found to be significant but we initially expect to find given this is a cancer dataset: P53, NOTCH_SIGNALING, DNA_REPAIR, HALLMARK_KRAS_SIGNALING_UP. First, we do not think there is an issue with the GSEA implementation, given what we reported above. Second, the original publication did not explicitly mention the p53 pathway. Last, [Tang et. al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6036325/) found the p53 pathway to be significantly enriched using 80 tumor and 20 normal samples. And [Gibbons et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3925633/) found that p53 mutations was high (46%) in non-small cell lung adenocarcinoma patients. We suspect the reason for the p53 pathway not being significant in this case might be due to the small sample size here (6 paired samples). So for now we will move forward with our GSEA analysis.

# In[ ]:


# Create "<local_dir>/GSEA_stats/" subdirectory
os.makedirs(os.path.join(local_dir, "GSEA_stats"), exist_ok=True)


# In[ ]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_filename -i num_runs -i statistic', '\nsource(\'../generic_expression_patterns_modules/GSEA_analysis.R\')\n\n# New files created: "<local_dir>/GSEA_stats/GSEA_stats_simulated_data_<project_id>_<n>.txt"\nfor (i in 0:(num_runs-1)) {\n    simulated_DE_stats_file <- paste(local_dir, \n                                     "DE_stats/DE_stats_simulated_data_", \n                                     project_id,\n                                     "_", \n                                     i,\n                                     ".txt",\n                                     sep = "")\n    \n    out_file <- paste(local_dir, \n                     "GSEA_stats/GSEA_stats_simulated_data_",\n                     project_id,\n                     "_",\n                     i,\n                     ".txt", \n                     sep = "")\n    \n    enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_filename, statistic) \n    \n    # Remove column with leading edge since its causing parsing issues\n    write.table(as.data.frame(enriched_pathways[1:7]), file = out_file, row.names = F, sep = "\\t")\n}')


# ### Rank pathways 

# In[ ]:


# Concatenate simulated experiments
simulated_GSEA_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id, 'GSEA')
simulated_GSEA_stats_all.set_index('pathway', inplace=True)
print(simulated_GSEA_stats_all.shape)


# In[ ]:


# Aggregate statistics across all simulated experiments
simulated_GSEA_summary_stats = calc.aggregate_stats(
    col_to_rank_pathways,
    simulated_GSEA_stats_all,
    'GSEA'
)

simulated_GSEA_summary_stats.head()


# In[ ]:


# Load association statistics for template experiment
template_GSEA_stats = template_enriched_pathways.iloc[:, :-1]
template_GSEA_stats.set_index('pathway', inplace=True)

template_GSEA_stats.head()

# Rank genes in template experiment
template_GSEA_stats = calc.rank_genes_or_pathways(
    col_to_rank_pathways,
    template_GSEA_stats,
    True
)


# In[ ]:


# Rank genes in simulated experiments
simulated_GSEA_summary_stats = calc.rank_genes_or_pathways(
    col_to_rank_pathways,
    simulated_GSEA_summary_stats,
    False
)


# ### Pathway summary table

# In[ ]:


# Create intermediate file: "<local_dir>/gene_summary_table_<col_to_rank_pathways>.tsv"
summary_pathway_ranks = process.generate_summary_table(
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
    col_to_rank_pathways,
    local_dir
)

summary_pathway_ranks.head()


# In[ ]:


# Create `pathway_summary_filename`
summary_pathway_ranks.to_csv(pathway_summary_filename, sep='\t')


# ### Compare gene ranking
# Studies have found that some genes are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represent a more systematic change. 
# 
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# In[ ]:


if compare_genes:
    # Get generic genes identified by Crow et. al.
    DE_prior_file = params['reference_gene_filename']
    ref_gene_col = params['reference_gene_name_col']
    ref_rank_col = params['reference_rank_col']
    
    figure_filename = os.path.join(
        local_dir, 
        f"gene_ranking_{col_to_rank_genes}.svg"
    )
    
    process.compare_gene_ranking(
        summary_gene_ranks,
        DE_prior_file,
        ref_gene_col,
        ref_rank_col,
        figure_filename
    )


# **Takeaway:**
# Based on the correlation plot, we can see that our simulation method is very good at capturing variability in genes that are very low or very high in the DE rank (i.e. are significantly differentially expressed often across different studies). These results serve to validate that our method can be used to identify these generic genes, as we were able to recapitulate some of the generic genes as those identified by Crow et. al. Additionally, our method extends the Crow et. al. work, which used array data, and since here we used RNA-seq.

# ### Compare pathway ranking try 1

# Studies have found that there are some pathways (gene sets) that are more likely to be significantly enriched in DEGs across a wide range of experimental designs. These generic pathways are not necessarily specific to the biological process being studied but instead represents a more systematic change.
# 
# We want to compare the ability to detect these generic pathways using our method vs those found by [Powers et. al.](https://www.biorxiv.org/content/10.1101/259440v1.full.pdf) publication.  We will use the `Hallmarks_qvalues_GSEAPreranked.csv` file from https://www.synapse.org/#!Synapse:syn11806255 as a reference. The file contains the q-value (adjusted p-value) for the test: given the enrichment score (ES) of the experiment is significant compared to the null distribution of enrichment scores, where the null set is generated from permuted gene sets. For each gene set (pathway) they calculate the q-value using this test. 
# 
# 
# To get a `reference ranking`, we calculate the fraction of experiments that a given pathway was significant (q-value <0.05) and use this rank pathways. `Our ranking` is to rank pathways based on the median q-value across the simulated experiments. We can then compare `our ranking` versus the `reference ranking.`

# In[ ]:


# Load Powers et. al. results file
powers_rank_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked.csv"
)


# In[ ]:


# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_filename, header=0, index_col=0)
powers_rank_df.drop(['Category'], axis=1, inplace=True)
print(powers_rank_df.shape)
powers_rank_df.head()


# In[ ]:


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
powers_rank_stats_df.head()


# In[ ]:


# Save reference file for input into comparison
powers_rank_processed_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked_processed.tsv"
)

powers_rank_stats_df.to_csv(powers_rank_processed_filename, sep="\t", )


# In[ ]:


if compare_genes:
    process.compare_pathway_ranking(
        summary_pathway_ranks,
        powers_rank_processed_filename,
    )


# The above shows that there is no correlation between our ranking (where pathways were ranked using median adjusted p-value score across simulated experiments) vs Powers et. al. ranking (where pathways were ranked based on the fraction of experiments they had adjusted p-value < 0.05). This is using the same workflow used to compare ranking of genes. Next let's try to use the fraction of adjusted p-value < 0.05 for our method and re-compare.

# In[ ]:


# Rank pathways based on the proportion of times they appear significant
simulated_GSEA_stats_all["significant"] = simulated_GSEA_stats_all['padj'] < 0.05
simulated_GSEA_stats_all_new = (
    simulated_GSEA_stats_all.groupby(
        ['pathway']
    )["significant"].sum() / num_runs
).to_frame(name="Fraction enriched")
simulated_GSEA_stats_all_new["Rank (simulated)"] = simulated_GSEA_stats_all_new.rank()

print(simulated_GSEA_stats_all_new.shape)
simulated_GSEA_stats_all_new.head()


# In[ ]:


if compare_genes:
    process.compare_pathway_ranking(
        simulated_GSEA_stats_all_new,
        powers_rank_processed_filename,
    )


# **Conclusion:**
# * If we compare the our ranking (`Rank (simulated)` column of the `summary_pathway_ranks` dataframe, we see that our highlight ranked pathways (rank >30) are consistent with those found to be generic (pathways found to be significantly enriched in >20% of experiments, figure 4A) in [Powers et. al.](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793) (EF2, TNFA, MYC_TARGET_V1/2, P53, HYPOXIA, INFLAMMATORY, APOPTOSIS, COAGULATION, KRAS) 
# 
# * Despite the eye ball consistency above, there is not a correlation between our method ranking and Powers et. al. ranking. The comparison we're doing here is not a precise match because our ranking is ES(pathway) from a null set while the Powers et. al. ranking is based on ES(pathway) vs null set. So the Powers et. al. ranking is corrected for my this null set. Though the comparision is not ideal we'd still expect a correlation in the ranking
# 
# Instead we will try to use the normalized ES (NES) values to rank the pathways for this next try.

# ## Compare rank pathways try 2

# In[ ]:


# Take absolute value of logFC and t statistic
simulated_GSEA_stats_all = process.abs_value_stats(simulated_GSEA_stats_all)


# In[ ]:


# Aggregate statistics across all simulated experiments
simulated_GSEA_summary_stats = calc.aggregate_stats(
    'NES',
    simulated_GSEA_stats_all,
    'GSEA'
)

simulated_GSEA_summary_stats.head()


# In[ ]:


# Take absolute value of NES
template_GSEA_stats = process.abs_value_stats(template_GSEA_stats)

# Rank genes in template experiment
template_GSEA_stats = calc.rank_genes_or_pathways(
    'NES',
    template_GSEA_stats,
    True
)


# In[ ]:


# Rank genes in simulated experiments
simulated_GSEA_summary_stats = calc.rank_genes_or_pathways(
    'NES',
    simulated_GSEA_summary_stats,
    False
)


# In[ ]:


summary_pathway_ranks = process.generate_summary_table(
    template_GSEA_stats,
    simulated_GSEA_summary_stats,
    'NES',
    local_dir
)

summary_pathway_ranks.head()


# In[ ]:


# Load Powers et. al. results file
powers_rank_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_NES_GSEAPreranked.csv"
)


# In[ ]:


# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_filename, header=0, index_col=0)
powers_rank_df.drop(['Category'], axis=1, inplace=True)
print(powers_rank_df.shape)
powers_rank_df.head()


# In[ ]:


# Rank pathways by NES score per experiment
# Get median rank per pathway
powers_rank_stats_df = powers_rank_df.abs().rank().median(axis=1).to_frame('Powers Rank')

powers_rank_stats_df.head()


# In[ ]:


# Save reference file for input into comparison
powers_rank_processed_filename = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_NES_GSEAPreranked_processed.tsv"
)

powers_rank_stats_df.to_csv(
    powers_rank_processed_filename, 
    sep="\t"
)


# In[ ]:


if compare_genes:
    process.compare_pathway_ranking(
        summary_pathway_ranks,
        powers_rank_processed_filename,
    )


# Try ranking pathways by NES for each simulated experiment and then having the final rank be the median(rank across simulated experiments)

# In[ ]:


# Concatenate simulated experiments
simulated_GSEA_stats_all = process.concat_simulated_data_columns(
    local_dir, 
    num_runs, 
    project_id, 
    'GSEA'
)
print(simulated_GSEA_stats_all.shape)
simulated_GSEA_stats_all.head()


# In[ ]:


# Rank pathways by NES score per experiment
# Get median rank per pathway
simulated_rank_stats_df = simulated_GSEA_stats_all.abs().rank().median(axis=1).to_frame('Rank (simulated)')

simulated_rank_stats_df.head()


# In[ ]:


if compare_genes:
    process.compare_pathway_ranking(
        simulated_rank_stats_df,
        powers_rank_processed_filename,
    )


# **Conclusion:**
# * Our top 20 pathways (out of a total of 50 pathways) are consistent with those found in Powers et. al. to be generic. If we don't consider the ranking of the pathways, our top 20 pathways overlap with those found to be frequently enriched in figure 4A of the Powers et. al. publication.
# * However, there is no correlation between our ranking vs Powers et. al ranking using adjusted p-values or NES values
# * This lack of correlation doesn't make sense to me, but not sure other things to try
# 
# Already verified:
# * Confirming that high rank = common pathways = low adjusted p-value, high NES
# * Confirm that high rank in our method is consistent with high rank in Powers et. al.
# * Increased number of simulated experiments to 100 
# * Verify that absolute value was taken before stats aggregated across simulated
