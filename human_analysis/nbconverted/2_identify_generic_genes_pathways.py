
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
col_to_rank = params['col_to_rank']
compare_genes = params['compare_genes']
statistic = params['gsea_statistic']

gene_summary_file = os.path.join(
    base_dir, 
    dataset_name, 
    f"generic_gene_summary_{project_id}.tsv")

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


# ### Simulate experiments using selected template experiment

# In[4]:


# Simulate multiple experiments
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
        i)


# Since this experiment contains both RNA-seq and smRNA-seq samples which are in different ranges so we will drop smRNA samples so that samples are within the same range. The analysis identifying these two subsets of samples can be found in this [notebook](../explore_data/0_explore_input_data.ipynb)

# In[5]:


if os.path.exists(sample_id_metadata_file):
    # Read in metadata
    metadata = pd.read_csv(sample_id_metadata_file, sep='\t', header=0, index_col=0)
    
    # Get samples to be dropped
    sample_ids_to_drop = list(metadata[metadata["processing"] == "drop"].index)

    process.subset_samples(sample_ids_to_drop,
                           num_runs,
                           local_dir,
                           project_id)


# In[6]:


# Round compendium read counts to int
#process.recast_int(num_runs, local_dir, project_id)


# ### Differential expression analysis

# In[7]:


# Load metadata file with grouping assignments for samples
metadata_file = os.path.join(
    base_dir,
    dataset_name,
    "data",
    "metadata",
    project_id+"_groups.tsv")


# In[8]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("DESeq2")')


# In[9]:


get_ipython().run_cell_magic('R', '', '# Load the DESeq2 library\nsuppressPackageStartupMessages(library("DESeq2"))')


# In[10]:


# Check ordering of sample ids is consistent between gene expression data and metadata
process.check_sample_ordering(template_data_file, metadata_file)


# In[11]:


get_ipython().run_cell_magic('R', '-i metadata_file -i project_id -i template_data_file -i local_dir', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\nget_DE_stats_DESeq(metadata_file,\n             project_id, \n             template_data_file,\n             "template",\n             local_dir,\n             "real")')


# In[12]:


# Check number of DEGs
template_DE_stats_file = os.path.join(
        local_dir,
        "DE_stats",
        f"DE_stats_template_data_{project_id}_real.txt")

template_DE_stats = pd.read_csv(template_DE_stats_file, sep="\t", header=0, index_col=0)

selected = template_DE_stats[(template_DE_stats['padj']<0.01) & (abs(template_DE_stats['log2FoldChange'])>1)]
print(selected.shape)


# In[13]:


# Check ordering of sample ids is consistent between gene expression data and metadata
for i in range(num_runs):
    simulated_data_file = os.path.join(
        local_dir,
        "pseudo_experiment",
        f"selected_simulated_data_{project_id}_{i}.txt")
        
    process.check_sample_ordering(simulated_data_file, metadata_file)


# In[14]:


get_ipython().run_cell_magic('R', '-i metadata_file -i project_id -i base_dir -i local_dir -i num_runs', '\nsource(\'../generic_expression_patterns_modules/DE_analysis.R\')\n\nnum_sign_DEGs_simulated <- c()\n\nfor (i in 0:(num_runs-1)){\n    simulated_data_file <- paste(local_dir, \n                                 "pseudo_experiment/selected_simulated_data_",\n                                 project_id,\n                                 "_", \n                                 i,\n                                 ".txt",\n                                 sep="")\n    \n    get_DE_stats_DESeq(metadata_file,\n                       project_id, \n                       simulated_data_file,\n                       "simulated",\n                       local_dir,\n                       i)\n}')


# **Validation:**
# * As a quick validation, [Kim et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3566005/) found 1459 DEGs (543 upregulated and 916 downregulated) using used the Bowtie and NEUMA applications for the mapping and quantification of RNA-Seq data. They used *edgeR* with a rigorous filtering procedure based on false discovery rates, minimum applicable patient numbers, and gene expression levels was devised to select reliable sets of DEGs and DEIs (see File S8 for details). For the
# 
# * Our results found ~3K DEGs which is close enough in range given that the data was processed using different methods. recount2 resource were aligned with the splice-aware Rail-RNA aligner

# ### Rank genes

# In[15]:


# Concatenate simulated experiments
simulated_DE_stats_all = process.concat_simulated_data(local_dir, num_runs, project_id)

print(simulated_DE_stats_all.shape)


# In[16]:


# Take absolute value of logFC and t statistic
simulated_DE_stats_all = process.abs_value_stats(simulated_DE_stats_all)


# In[17]:


# Aggregate statistics across all simulated experiments
simulated_DE_summary_stats = calc.aggregate_stats(col_to_rank,
                                                  simulated_DE_stats_all)


# In[18]:


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
template_DE_stats = calc.rank_genes(col_to_rank,
                                   template_DE_stats,
                                   True)


# In[19]:


# Rank genes in simulated experiments
simulated_DE_summary_stats = calc.rank_genes(col_to_rank,
                                            simulated_DE_summary_stats,
                                            False)


# ### Gene summary table

# In[20]:


summary_gene_ranks = process.generate_summary_table(template_DE_stats,
                                                   simulated_DE_summary_stats,
                                                   col_to_rank,
                                                   local_dir)

summary_gene_ranks.head()


# In[21]:


summary_gene_ranks.to_csv(
    gene_summary_file, sep='\t')


# ### GSEA 
# **Goal:** To detect modest but coordinated changes in prespecified sets of related genes (i.e. those genes in the same pathway or share the same GO term).
# 
# 1. Ranks all genes based using DE association statistics. In this case we used the p-value scores to rank genes. logFC returned error -- need to look into this.
# 2. An enrichment score (ES) is defined as the maximum distance from the middle of the ranked list. Thus, the enrichment score indicates whether the genes contained in a gene set are clustered towards the beginning or the end of the ranked list (indicating a correlation with change in expression). 
# 3. Estimate the statistical significance of the ES by a phenotypic-based permutation test in order to produce a null distribution for the ES( i.e. scores based on permuted phenotype)

# In[22]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("GSA")\n#BiocManager::install("fgsea")')


# In[23]:


get_ipython().run_cell_magic('R', '', 'suppressPackageStartupMessages(library("GSA"))\nsuppressPackageStartupMessages(library("fgsea"))')


# In[24]:


# Load pathway data
hallmark_DB_file = os.path.join(
    local_dir,
    "hallmark_DB.gmt")


# In[25]:


get_ipython().run_cell_magic('R', '-i template_DE_stats_file -i hallmark_DB_file -i statistic', '# This cell is temporary and contains the same code as the GSEA_analysis.R file\n# This is being used to debug the GSEA results found\n# Read in data\nDE_stats_data <- read.table(template_DE_stats_file, sep="\\t", header=TRUE, row.names=NULL)\n\n# feature 1: numeric vector\nif (statistic == \'logFC\'){\n    col_num = 2\n} else if (statistic == \'log2FoldChange\'){\n    col_num = 3\n} else if (statistic ==\'t\'){\n    col_num = 4\n} else if (statistic == \'p-value\'){\n    col_num = 5\n} else if (statistic == \'adj p-value\' || statistic == \'pvalue\'){\n    col_num = 6\n} else if ( statistic == \'padj\'){\n    col_num = 7\n}\n\nrank_genes <- as.numeric(as.character(DE_stats_data[,col_num]))\n\n# feature 2: named vector of gene ids\nnames(rank_genes) <- as.character(DE_stats_data[,1])\nprint(head(rank_genes))\n## feature 3: decreasing order\nrank_genes <- sort(rank_genes, decreasing = TRUE)\nprint(head(rank_genes))\n\npathway_DB_data <- gmtPathways(hallmark_DB_file)\nprint(head(pathway_DB_data))\n\n#print(head(pathway_DB_data))\n# GSEA is a generic gene set enrichment function\n# Different backend methods can be applied depending on the \n# type of annotations\n# Here we will use fgsea\n\n#enrich_pathways <- fgsea(pathways=pathway_DB_data,\n#                         stats=rank_genes,\n#                         nperm=10000,\n#                         minSize=10,\n#                         maxSize=200\n#                        )\n\n#print(enrich_pathways)\n\n# Look at those pathways WITH significant padj\n#plotEnrichment(pathway_DB_data[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)\n#plotEnrichment(pathway_DB_data[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)\n\n# Look at those pathways with NOT significant padj\nplotEnrichment(pathway_DB_data[["HALLMARK_P53_PATHWAY"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)\n#barplot(sort(rank_genes, decreasing = T))"""')


# In[26]:


get_ipython().run_cell_magic('R', '-i template_DE_stats_file -i hallmark_DB_file -i statistic -o template_enriched_pathways', "\nsource('../generic_expression_patterns_modules/GSEA_analysis.R')\ntemplate_enriched_pathways <- find_enriched_pathways(template_DE_stats_file, hallmark_DB_file, statistic)")


# ## Quick validation
# Significant pathways found include: HALLMARK_G2M_CHECKPOINT, TNFA_SIGNALING, INFLAMMATORY_RESPONSE, APOPTOSIS, HYPOXIA are expected and consistent with what the publication found. However, there are pathways not found to be significant that we would expect to find: P53, NOTCH_SIGNALING, DNA_REPAIR. Not sure how to explain this since the enrichment score plot looks consistent with the p-value score .

# In[27]:


print(template_enriched_pathways.shape)
template_enriched_pathways[template_enriched_pathways['padj'] < 0.05].sort_values(by='padj')


# In[28]:


template_enriched_pathways[template_enriched_pathways['padj'] >= 0.05].sort_values(by='padj')


# In[29]:


template_enriched_pathways.set_index('pathway', inplace=True)
template_enriched_pathways.loc["HALLMARK_P53_PATHWAY"]


# In[30]:


# Check expression of TP53 genes in cancer vs normal samples
tp53_genes = list(template_DE_stats[list(template_DE_stats.index.str.find("TP53") == 0)].index)
template_data = pd.read_csv(template_data_file, sep="\t", index_col=0)

selected_template_data = template_data[tp53_genes]


# In[31]:


# Load metadata file with grouping assignments for samples
metadata_file = os.path.join(
    "data",
    "metadata",
    project_id+"_groups.tsv")

# Read metadata
metadata = pd.read_csv(
    metadata_file,
    header=0,
    sep='\t',
    index_col=0)


# In[32]:


# PCA encode
pca = PCA(n_components=2)

model = pca.fit(selected_template_data)
template_PCAencoded = model.transform(selected_template_data)

template_PCAencoded_df = pd.DataFrame(data=template_PCAencoded ,
                                         index=template_data.index,
                                         columns=['1','2'])


# In[33]:


# Add tumor/normal labels
template_data_labeled = pd.merge(template_PCAencoded_df, metadata, left_index=True, right_index=True)
template_data_labeled.loc[template_data_labeled['group'] == 1,'source'] = 'Normal'
template_data_labeled.loc[template_data_labeled['group'] == 2,'source'] = 'Tumor'
template_data_labeled.head()


# In[34]:


# Plot
fig = ggplot(template_data_labeled, aes(x='1', y='2'))
fig += geom_point(aes(color='source'), alpha=0.7)
fig += labs(x ='PCA 1',
            y = 'PCA 2',
            title = 'PCA template data')
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


# In[37]:


get_ipython().run_cell_magic('R', '-i project_id -i local_dir -i hallmark_DB_file -i num_runs -i statistic', '\nsource(\'../generic_expression_patterns_modules/GSEA_analysis.R\')\n\nfor (i in 0:(num_runs-1)){\n    simulated_DE_stats_file <- paste(local_dir, \n                                 "DE_stats/DE_stats_simulated_data_", \n                                 project_id,\n                                 "_", \n                                 i,\n                                 ".txt",\n                                 sep="")\n    \n    out_file = paste(local_dir, \n                     "GSEA_stats/GSEA_simulated_data_",\n                     project_id,\n                     "_",\n                     i,\n                     ".txt", \n                     sep="")\n    \n    enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_file, statistic) \n    \n    write.table(data.frame(enriched_pathways), file = out_file, row.names = T, sep = "\\t")\n    }')


# ### Rank pathways 

# ### Pathway summary table

# ### Compare gene ranking
# Studies have found that there are some genes that are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represents a more systematic change. 
# 
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# In[36]:


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
        
    # Drop genes with 0 mean base expression
    # Note: These lowly expressed genes were not pre-filtered before DESeq
    # (Michael Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of dispersion is difficult 
    # for genes with very, very low average counts. See the methods. 
    # However it doesn't really effect the outcome because these genes have almost no power for detecting 
    # differential expression. Effects runtime though.
    shared_gene_rank_scaled_df = shared_gene_rank_scaled_df[~shared_gene_rank_scaled_df['Rank (simulated)'].isna()]
    
    # Get correlation
    r, p, ci_high, ci_low = calc.spearman_ci(0.95,
                                             shared_gene_rank_scaled_df,
                                             1000)
    print(r, p, ci_high, ci_low)
    
    # Plot our ranking vs published ranking
    fig_file = os.path.join(
        local_dir, 
        "gene_ranking_"+col_to_rank+".svg")

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
