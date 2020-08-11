
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
process.recast_int(num_runs, local_dir, project_id)


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


# Based on publication, they found most enriched terms include cancer, cellular growth, cell proliferation, cell death
# 
# Would would expect similar enrichment in these types of pathways

# In[25]:


"""%%R -i template_DE_stats_file -i hallmark_DB_file -i statistic
# Read in data
DE_stats_data <- read.table(template_DE_stats_file, sep="\t", header=TRUE, row.names=NULL)

# Sort genes by feature 1
statistic = 't'

# feature 1: numeric vector
if (statistic =='t'){
    col_num = 4
} else if (statistic == 'adj p-value'){
    col_num = 6
} else if (statistic == 'p-value'){
    col_num = 5
} else if (statistic == 'logFC'){
    col_num = 2
}

rank_genes <- as.numeric(as.character(DE_stats_data[,col_num]))

# feature 2: named vector of gene ids
names(rank_genes) <- as.character(DE_stats_data[,1])

## feature 3: decreasing order
rank_genes <- sort(rank_genes, decreasing = TRUE)
print(head(rank_genes))

pathway_DB_data <- gmtPathways(hallmark_DB_file)
print(head(pathway_DB_data))
#pathway_parsed <- {}
#for (i in 1:length(pathway_DB_data$genesets)){
#pathway_parsed[pathway_DB_data$geneset.name[i]] <- as.list(pathway_DB_data$genesets[i])
#}

#print(head(pathway_DB_data))
# GSEA is a generic gene set enrichment function
# Different backend methods can be applied depending on the 
# type of annotations
# Here we will use fgsea
#enrich_pathways <- GSEA(geneList=rank_genes, 
#                        TERM2GENE=pathway_DB_data,
#                        nPerm=100000,
#                        by='fgsea',
#                        verbose=T)
enrich_pathways <- fgsea(pathways=pathway_DB_data,
                         stats=rank_genes,
                         nperm=10000,
                         minSize=10,
                         maxSize=200
                        )

print(enrich_pathways)
#plotEnrichment(pathway_parsed[["HALLMARK_DNA_REPAIR"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)
#plotEnrichment(pathway_parsed[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)
#plotEnrichment(pathway_parsed[["HALLMARK_P53_PATHWAY"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)
#plotEnrichment(pathway_DB_data[["HALLMARK_BILE_ACID_METABOLISM"]], stats=rank_genes, gseaParam = 1, ticksSize = 0.2)
#barplot(sort(rank_genes, decreasing = T))"""


# In[26]:


get_ipython().run_cell_magic('R', '-i template_DE_stats_file -i hallmark_DB_file -i statistic -o template_enriched_pathways', "\nsource('../generic_expression_patterns_modules/GSEA_analysis.R')\ntemplate_enriched_pathways <- find_enriched_pathways(template_DE_stats_file, hallmark_DB_file, statistic)")


# In[29]:


print(template_enriched_pathways.shape)
#template_enriched_pathways['size'].max()
template_enriched_pathways[template_enriched_pathways['padj'] < 0.05]
#template_enriched_pathways[template_enriched_pathways['padj']<0.1]
#template_enriched_pathways.head()


# In[30]:


template_enriched_pathways.set_index('pathway', inplace=True)
template_enriched_pathways.loc["HALLMARK_P53_PATHWAY"]


# In[ ]:


"""%%R -i project_id -i local_dir -i hallmark_DB_file -i num_runs -i statistic

source('../generic_expression_patterns_modules/GSEA_analysis.R')

for (i in 0:(num_runs-1)){
    simulated_DE_stats_file <- paste(local_dir, 
                                 "DE_stats/DE_stats_simulated_data_", 
                                 project_id,
                                 "_", 
                                 i,
                                 ".txt",
                                 sep="")
    
    out_file = paste(local_dir, 
                     "GSEA_stats/GSEA_simulated_data_",
                     project_id,
                     "_",
                     i,
                     ".txt", 
                     sep="")
    
    enriched_pathways <- find_enriched_pathways(simulated_DE_stats_file, hallmark_DB_file, statistic) 
    #print(head(enriched_pathways))
    
    write.table(enriched_pathways, file = out_file, row.names = T, sep = "\t")
    }"""


# ### Rank pathways 

# ### Pathway summary table

# ### Compare gene ranking
# Studies have found that there are some genes that are more likely to be differentially expressed even across a wide range of experimental designs. These *generic genes* are not necessarily specific to the biological process being studied but instead represents a more systematic change. 
# 
# We want to compare the ability to detect these generic genes using our method vs those found by [Crow et. al. publication](https://www.pnas.org/content/pnas/116/13/6491.full.pdf). Their genes are ranked 0 = not commonly DE; 1 = commonly DE. Genes by the number differentially expressed gene sets they appear in and then ranking genes by this score.

# In[48]:


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
    # (Micheal Love, author of DESeq2): In our DESeq2 paper we discuss a case where estimation of dispersion is difficult 
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
