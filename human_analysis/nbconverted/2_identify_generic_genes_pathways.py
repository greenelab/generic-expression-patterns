
# coding: utf-8

# # Identify generic genes and pathways
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
import sys
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import pickle
from sklearn.decomposition import PCA

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
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import calc, process

np.random.seed(123)


# In[2]:


# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

config_file = os.path.abspath(os.path.join(base_dir,
                                           "configs",
                                           "config_human.tsv"))
params = utils.read_config(config_file)


# In[3]:


# Load params
local_dir = params["local_dir"]
dataset_name = params['dataset_name']
NN_architecture = params['NN_architecture']
num_runs = params['num_simulated']
project_id = params['project_id']
metadata_col_id = params['metadata_colname']
template_data_file = params['template_data_file']
original_compendium_file = params['compendium_data_file']
normalized_compendium_file = params['normalized_compendium_data_file']
scaler_file = params['scaler_transform_file']
col_to_rank_genes = params['rank_genes_by']
col_to_rank_pathways = params['rank_pathways_by']
compare_genes = params['compare_genes']
statistic = params['gsea_statistic']

NN_dir = os.path.join(
    base_dir, 
    dataset_name, 
    "models", 
    NN_architecture)

# Load metadata file with grouping assignments for samples
sample_id_metadata_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    f"{project_id}_process_samples.tsv")

# Load pickled file
scaler = pickle.load(open(scaler_file, "rb"))


# In[4]:


# Output files
gene_summary_file = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_gene_summary_{project_id}.tsv")

pathway_summary_file = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_pathway_summary_{project_id}.tsv")


# ### Simulate experiments using selected template experiment

# In[5]:


"""# Simulate multiple experiments
for i in range(num_runs):
    simulate_expression_data.shift_template_experiment(
        normalized_compendium_file,
        project_id,
        metadata_col_id,
        NN_architecture,
        dataset_name,
        scaler,
        local_dir,
        base_dir,
        i)"""


# Since this experiment contains both RNA-seq and smRNA-seq samples which are in different ranges so we will drop smRNA samples so that samples are within the same range. The analysis identifying these two subsets of samples can be found in this [notebook](../explore_data/0_explore_input_data.ipynb)

# In[6]:


"""if os.path.exists(sample_id_metadata_file):
    # Read in metadata
    metadata = pd.read_csv(sample_id_metadata_file, sep='\t', header=0, index_col=0)
    
    # Get samples to be dropped
    sample_ids_to_drop = list(metadata[metadata["processing"] == "drop"].index)

    process.subset_samples(sample_ids_to_drop,
                           num_runs,
                           local_dir,
                           project_id)"""


# In[7]:


# Round simulated read counts to int in order to run DESeq
process.recast_int(num_runs, local_dir, project_id)


# ### Differential expression analysis

# In[8]:


# Load metadata file with grouping assignments for samples
metadata_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    project_id+"_groups.tsv")


# In[9]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("DESeq2")')


# In[10]:


get_ipython().run_cell_magic('R', '', '# Load the DESeq2 library\nsuppressPackageStartupMessages(library("DESeq2"))')


# In[11]:


# Check ordering of sample ids is consistent between gene expression data and metadata
process.compare_and_reorder_samples(template_data_file, metadata_file)


# In[12]:


get_ipython().run_cell_magic('R', '-i metadata_file -i project_id -i template_data_file -i local_dir', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\nget_DE_stats_DESeq(metadata_file,\n             project_id, \n             template_data_file,\n             "template",\n             local_dir,\n             "real")')


# In[13]:


# Check number of DEGs
template_DE_stats_file = os.path.join(
        local_dir,
        "DE_stats",
        f"DE_stats_template_data_{project_id}_real.txt")

template_DE_stats = pd.read_csv(template_DE_stats_file, sep="\t", header=0, index_col=0)

selected = template_DE_stats[(template_DE_stats['padj']<0.01) & (abs(template_DE_stats['log2FoldChange'])>1)]
print(selected.shape)


# In[14]:


# Check ordering of sample ids is consistent between gene expression data and metadata
for i in range(num_runs):
    simulated_data_file = os.path.join(
        local_dir,
        "pseudo_experiment",
        f"selected_simulated_data_{project_id}_{i}.txt")
        
    process.compare_and_reorder_samples(simulated_data_file, metadata_file)


# In[15]:


get_ipython().run_cell_magic('R', '-i metadata_file -i project_id -i base_dir -i local_dir -i num_runs', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\nnum_sign_DEGs_simulated <- c()\n\nfor (i in 0:(num_runs-1)){\n    simulated_data_file <- paste(local_dir, \n                                 "pseudo_experiment/selected_simulated_data_",\n                                 project_id,\n                                 "_", \n                                 i,\n                                 ".txt",\n                                 sep="")\n    \n    get_DE_stats_DESeq(metadata_file,\n                       project_id, \n                       simulated_data_file,\n                       "simulated",\n                       local_dir,\n                       i)\n}')


# **Validation:**
# * As a quick validation, [Kim et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3566005/) found 1459 DEGs (543 upregulated and 916 downregulated) using used the Bowtie and NEUMA applications for the mapping and quantification of RNA-Seq data. They used *edgeR* with a rigorous filtering procedure based on false discovery rates, minimum applicable patient numbers, and gene expression levels was devised to select reliable sets of DEGs and DEIs (see File S8 for details). For the
# 
# * Our results found ~3K DEGs which is close enough in range given that the data was processed using different methods. recount2 resource were aligned with the splice-aware Rail-RNA aligner

# ### Rank genes

# In[16]:


# Concatenate simulated experiments
simulated_DE_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id, 'DE')

print(simulated_DE_stats_all.shape)


# In[17]:


# Take absolute value of logFC and t statistic
simulated_DE_stats_all = process.abs_value_stats(simulated_DE_stats_all)


# In[18]:


# Aggregate statistics across all simulated experiments
simulated_DE_summary_stats = calc.aggregate_stats(col_to_rank_genes,
                                                  simulated_DE_stats_all,
                                                  'DE')


# In[19]:


# Load association statistics for template experiment
template_DE_stats_file = os.path.join(
    local_dir,
    "DE_stats",
    "DE_stats_template_data_"+project_id+"_real.txt")

template_DE_stats = pd.read_csv(
    template_DE_stats_file,
    header=0,
    sep='\t',
    index_col=0)

# Take absolute value of logFC and t statistic
template_DE_stats = process.abs_value_stats(template_DE_stats)

# Rank genes in template experiment
template_DE_stats = calc.rank_genes_or_pathways(col_to_rank_genes,      
                                                template_DE_stats,
                                                True)


# In[20]:


# Rank genes in simulated experiments
simulated_DE_summary_stats = calc.rank_genes_or_pathways(col_to_rank_genes,
                                                         simulated_DE_summary_stats,
                                                         False)


# ### Gene summary table

# In[21]:


summary_gene_ranks = process.generate_summary_table(template_DE_stats,
                                                   simulated_DE_summary_stats,
                                                   col_to_rank_genes,
                                                   local_dir)

summary_gene_ranks.head()


# In[22]:


summary_gene_ranks.to_csv(
    gene_summary_file, sep='\t')


# ### GSEA 
# **Goal:** To detect modest but coordinated changes in prespecified sets of related genes (i.e. those genes in the same pathway or share the same GO term).
# 
# 1. Ranks all genes based using DE association statistics. In this case we used the p-value scores to rank genes. logFC returned error -- need to look into this.
# 2. An enrichment score (ES) is defined as the maximum distance from the middle of the ranked list. Thus, the enrichment score indicates whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list (indicating a correlation with change in expression). 
# 3. Estimate the statistical significance of the ES by a phenotypic-based permutation test in order to produce a null distribution for the ES( i.e. scores based on permuted phenotype)

# In[23]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("GSA")\n#BiocManager::install("fgsea")')


# In[24]:


get_ipython().run_cell_magic('R', '', 'suppressPackageStartupMessages(library("GSA"))\nsuppressPackageStartupMessages(library("fgsea"))')


# In[25]:


# Load pathway data
hallmark_DB_file = os.path.join(
    local_dir,
    "hallmark_DB.gmt")


# In[26]:


get_ipython().run_cell_magic('R', '-i template_DE_stats_file -i hallmark_DB_file -i statistic -o template_enriched_pathways', "\nsource('../generic_expression_patterns_modules/GSEA_analysis.R')\ntemplate_enriched_pathways <- find_enriched_pathways(template_DE_stats_file, hallmark_DB_file, statistic)")


# In[27]:


print(template_enriched_pathways.shape)
template_enriched_pathways[template_enriched_pathways['padj'] < 0.05].sort_values(by='padj').head()


# In[28]:


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

# In[29]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_file -i num_runs -i statistic', '\nsource(\'../generic_expression_patterns_modules/GSEA_analysis.R\')\n\nfor (i in 0:(num_runs-1)){\n    simulated_DE_stats_file <- paste(local_dir, \n                                 "DE_stats/DE_stats_simulated_data_", \n                                 project_id,\n                                 "_", \n                                 i,\n                                 ".txt",\n                                 sep="")\n    \n    out_file = paste(local_dir, \n                     "GSEA_stats/GSEA_stats_simulated_data_",\n                     project_id,\n                     "_",\n                     i,\n                     ".txt", \n                     sep="")\n    \n    enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_file, statistic) \n    \n    # Remove column with leading edge since its causing parsing issues\n    write.table(as.data.frame(enriched_pathways[1:7]), file = out_file, row.names = F, sep = "\\t")\n    }')


# ### Rank pathways 

# In[30]:


# Concatenate simulated experiments
simulated_GSEA_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id, 'GSEA')
simulated_GSEA_stats_all.set_index('pathway', inplace=True)
print(simulated_GSEA_stats_all.shape)


# In[31]:


simulated_GSEA_stats_all.head()


# In[32]:


simulated_GSEA_stats_all.groupby(["pathway"])


# In[33]:


# Aggregate statistics across all simulated experiments
simulated_GSEA_summary_stats = calc.aggregate_stats(col_to_rank_pathways,
                                                    simulated_GSEA_stats_all,
                                                    'GSEA')

simulated_GSEA_summary_stats.head()


# In[34]:


# Load association statistics for template experiment
template_GSEA_stats = template_enriched_pathways.iloc[:,:-1]
template_GSEA_stats.set_index('pathway', inplace=True)

template_GSEA_stats.head()

# Rank genes in template experiment
template_GSEA_stats = calc.rank_genes_or_pathways(col_to_rank_pathways,
                                                  template_GSEA_stats,
                                                  True)


# In[35]:


# Rank genes in simulated experiments
simulated_GSEA_summary_stats = calc.rank_genes_or_pathways(col_to_rank_pathways,
                                                           simulated_GSEA_summary_stats,
                                                           False)


# In[36]:


simulated_GSEA_summary_stats.head()


# In[37]:


template_GSEA_stats.head()


# ### Pathway summary table

# In[38]:


summary_pathway_ranks = process.generate_summary_table(template_GSEA_stats,
                                                       simulated_GSEA_summary_stats,
                                                       col_to_rank_pathways,
                                                       local_dir)

summary_pathway_ranks.head()


# In[39]:


summary_pathway_ranks.sort_values(by="Rank (simulated)", ascending=False)


# In[40]:


summary_pathway_ranks.to_csv(
    pathway_summary_file, sep='\t')


# ### Compare gene ranking
# Studies have found that there are some genes that are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represents a more systematic change. 
# 
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# In[41]:


if compare_genes:
    # Get generic genes identified by Crow et. al.
    DE_prior_file = params['reference_gene_file']
    ref_gene_col = params['reference_gene_name_col']
    ref_rank_col = params['reference_rank_col']
    
    # Merge our ranking and reference ranking
    shared_gene_rank_df = process.merge_ranks_to_compare(
        summary_gene_ranks,
        DE_prior_file,
        ref_gene_col,
        ref_rank_col)
    
    if max(shared_gene_rank_df["Rank (simulated)"]) != max(shared_gene_rank_df[ref_rank_col]):
        shared_gene_rank_scaled_df = process.scale_reference_ranking(shared_gene_rank_df, ref_rank_col)
        
    # Note: These lowly expressed genes were not pre-filtered before DESeq
    # (Micheal Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of dispersion is difficult 
    # for genes with very, very low average counts. See the methods. 
    # However it doesn't really effect the outcome because these genes have almost no power for detecting 
    # differential expression. Effects runtime though.
    shared_gene_rank_scaled_df = shared_gene_rank_scaled_df[~shared_gene_rank_scaled_df['Rank (simulated)'].isna()]
    
    # Get correlation
    r, p, ci_low, ci_high = calc.spearman_ci(0.95,
                                             shared_gene_rank_scaled_df,
                                             1000,
                                             'DE')
    print(r, p, ci_low, ci_high)
    
    # Plot our ranking vs published ranking
    fig_file = os.path.join(
        local_dir, 
        "gene_ranking_"+col_to_rank_genes+".svg")

    fig = sns.jointplot(data=shared_gene_rank_scaled_df,
                        x='Rank (simulated)',
                        y=ref_rank_col,
                        kind='hex',
                        marginal_kws={'color':'white'})
    fig.set_axis_labels("Our preliminary method", "DE prior (Crow et. al. 2019)", fontsize=14)

    fig.savefig(fig_file,
                format='svg',
                bbox_inches="tight",
                transparent=True,
                pad_inches=0,
                dpi=300,)


# **Takeaway:**
# Based on the correlation plot, we can see that our simulation method is very good at capturing variability in genes that are very low or very high in the DE rank (i.e. are significantly differentially expressed often across different studies). These results serve to validate that our method can be used to identify these generic genes, as we were able to recapitulate the some of the generic genes as those identified by Crow et. al. Additionally, our method extends the Crow et. al. work, which used array data, and since here we used RNA-seq.

# ### Compare pathway ranking

# Studies have found that there are some pathways (gene sets) that are more likely to be significantly enriched in DEGs across a wide range of experimental designs. These generic pathways are not necessarily specific to the biological process being studied but instead represents a more systematic change.
# 
# We want to compare the ability to detect these generic pathways using our method vs those found by [Powers et. al.](https://www.biorxiv.org/content/10.1101/259440v1.full.pdf) publication.  Significant gene sets are defined as an FDR < 0.05. Gene sets are ranked by the difference in
# the fraction of experiments with significant positive and negative enrichment
# 
# Using the `Hallmarks_qvalues_GSEAPreranked.csv` file from https://www.synapse.org/#!Synapse:syn11806255 as a reference.

# In[42]:


# Load Powers et. al. results file
powers_rank_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked.csv")


# In[43]:


# Read Powers et. al. data
# This file contains qvalue results for hallmark pathways across ~400 experiments
powers_rank_df = pd.read_csv(powers_rank_file, header=0, index_col=0)
powers_rank_df.drop(['Category'], axis=1, inplace=True)
print(powers_rank_df.shape)
powers_rank_df.head()


# In[44]:


# Count the number of experiments where a given pathway was found to be enriched (qvalue < 0.05)
total_num_experiments = powers_rank_df.shape[1]
frac_enriched_pathways = ((powers_rank_df< 0.05).sum(axis=1)/total_num_experiments)

# Rank pathways from 0-50, 50 indicating that the pathways was frequently enriched
pathway_ranks = frac_enriched_pathways.rank()

powers_rank_stats_df = pd.DataFrame(
    data={'Fraction enriched': frac_enriched_pathways.values,
          'Powers Rank':pathway_ranks.values
         },
    index=powers_rank_df.index
)
powers_rank_stats_df.head()


# In[45]:


# Save reference file for input into comparison
powers_rank_processed_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    "Hallmarks_qvalues_GSEAPreranked_processed.tsv")

powers_rank_stats_df.to_csv(powers_rank_processed_file, sep="\t", )


# In[46]:


if compare_genes:
    # Get column headers for generic pathways identified by Powers et. al.
    ref_gene_col='index'
    ref_rank_col = 'Powers Rank'
    
    # Merge our ranking and reference ranking
    shared_pathway_rank_df = process.merge_ranks_to_compare(
        summary_pathway_ranks,
        powers_rank_processed_file,
        ref_gene_col,
        ref_rank_col
        )
    
    if max(shared_pathway_rank_df["Rank (simulated)"]) != max(shared_pathway_rank_df[ref_rank_col]):
        shared_pathway_rank_scaled_df = process.scale_reference_ranking(shared_pathway_rank_df, ref_rank_col)
        
    # Note: These lowly expressed genes were not pre-filtered before DESeq
    # (Micheal Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of dispersion is difficult 
    # for genes with very, very low average counts. See the methods. 
    # However it doesn't really effect the outcome because these genes have almost no power for detecting 
    # differential expression. Effects runtime though.
    
    shared_pathway_rank_scaled_df = shared_pathway_rank_scaled_df[~shared_pathway_rank_scaled_df['Rank (simulated)'].isna()]
    
    # Get correlation
    r, p, ci_low, ci_high = calc.spearman_ci(0.95,
                                             shared_pathway_rank_scaled_df,
                                             1000,
                                             'GSEA')
    print(r, p, ci_low, ci_high)
    
    # Plot our ranking vs published ranking
    fig_file = os.path.join(
        local_dir, 
        "pathway_ranking_"+col_to_rank_pathways+".svg")

    fig = sns.jointplot(data=shared_pathway_rank_scaled_df,
                        x='Rank (simulated)',
                        y=ref_rank_col,
                        kind='hex',
                        marginal_kws={'color':'white'})
    fig.set_axis_labels("Our preliminary method", "Gene set rank (Powers et. al. 2018)", fontsize=14)

    fig.savefig(fig_file,
                format='svg',
                bbox_inches="tight",
                transparent=True,
                pad_inches=0,
                dpi=300)

