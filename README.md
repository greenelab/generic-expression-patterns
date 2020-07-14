# Rank pathways analysis

**Alexandra J Lee, James C Costello and Casey S Greene**

**May 2020**

**University of Pennsylvania, University of Colorado Anschutz Medical Campus**

**Background** (from [Powers et. al., Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793)) 
* Gene Set Enrichment Analysis (GSEA) was developed to help with the analysis and interpretation of the long lists of genes produced from high-throughput transcriptomic experiments. 
* By summarizing genome-wide gene expression changes into gene sets (groups of functionally related genes) a user can gain insight into how biological pathways and processes are affected under the tested experimental conditions. 
* The power of GSEA lies in its use of gene sets, which provide a more stable and interpretable measure of biological functions compared to individual genes that can show greater experimental and technical variation
* The underlying hypothesis of GSEA is that genes involved in a similar biological process or pathway (grouped into gene sets) are coordinately regulated. Thus, if an experimental perturbation activates a pathway, the genes in the associated gene set will be coordinately up-regulated (i.e. there will be an overrepresentation of genes in the gene set in the set of DEGs) and this pattern can be identified using statistical tests. The enrichment score, which reflects the degree to which genes in a gene set are over-represented at either end of a ranked gene list

**Rationale**: People performing differential expression (DE) analysis find that some genes and subsequent pathways are more likely to be differentially expressed even across a wide range of experimental designs.[Powers et. al., Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793 ); [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491). Powers et. al. developped a tool, [Explorer-InContext](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793) with corresponding [app](https://www.biorxiv.org/content/10.1101/659847v1.full.pdf), to try to correct for these commonly enriched gene sets by comparing gene set ranks from a target experiment with a null set of experiments (called "context"). In other words, the gene set ranks obtained from the target experiment are compared against the gene set ranks from the null experiments to determine if high rank gene sets from the target experiment are significant given the distribution of their rank in the null set. This method required a large manual curation effort to: collect a large set of samples with corresponding metadata (metadata is used to group samples per experiment and perform DE analysis to get ranked list of genes)

**Problem**: 
* This method required a large manual curation effort
* If you want to perform a new DE analysis in a different biological **context** (i.e. different organism, tissue, media) then you might not have the curated data available. Switch contexts will require a lot of manual effort. 
* Similarly, using a different statistical method will require re-curation effort

**Question**: Can we use our gene expression simulator, [ponyo](https://github.com/greenelab/ponyo), automatically generate null experiments for different contexts in order to overcome this manual effort?

## Initial experiment
**Question**: Can our compendia simulation method, [ponyo](https://github.com/greenelab/ponyo) identify the same commonly DEGs found in [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491) and commonly enriched pathways as filtered out in GSEA-InContext?

**Approach**:
1. Select treatment vs control experiment from recount2
2. Simulate 100 new experiments using experiment (1) as template
3. Perform DE analysis to get association statistics
4. Perform enrichment analysis using DE stats
5. Rank DEGs and pathways based on aggregated statistics across the simulated experiments

**Hypothesis**: 
* Our compendia simulation method, [ponyo](https://github.com/greenelab/ponyo) can help us identify the same commonly DEGs found in [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491) 
* If we ranked pathways by the number of times that they showed up in our ponyo-generated null set, and we looked at those, we'd find the pathways that [GSEA-InContext](https://www.biorxiv.org/content/10.1101/659847v1) is designed to not find.


## How to run notebooks from generic-expression-patterns

In order to run this simulation on your own gene expression data the following steps should be performed:

First you need to set up your local repository: 
1. Clone the `generic-expression-patterns` repository
2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
3. Set up conda environment by running the following command in the terminal:
```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate generic_expression

pip install -e .
```
4. Navigate to either the `human_analysis` or `pseudomonas_analysis` directories and run the notebooks.

## How to run using your own data

In order to run this simulation on your own gene expression data the following steps should be performed:

First you need to set up your local repository and environment: 
1. Clone the `generic-expression-patterns` repository
2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
3. Set up conda environment by running the following command in the terminal:
```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate generic_expression

pip install -e .
```
4. Create a new analysis folder in the main directory. This is equivalent to the `human_analysis` directory
5. Copy jupyter notebooks (1_process_data.ipynb, 2_identify_generic_genes_pathways.ipynb) into analysis directory.
6. Customize `1_process_data.ipynb` to generate the following saved files: 1) compendium of gene expression data, 2) template experiment data, 3) normalized gene expression compendium. See examples of this processing in `human_analysis/` and `pseudomonas_analysis/`.
7. Required data files:
* Gene expression compendium to use as input data. Specify path of data file in config file.
* Metadata matrix (sample x experimental metadata) with sample id as index. Add file to `analysis_dir/data/metadata/`. Specify path of metadata file in config file.
* Sample grouping matrix (sample x group) with group = [1 if control; 2 if case]. Add file to `analysis_dir/data/metadata/<project_id>_groups.tsv`.
8. `2_identify_generic_genes_pathways.ipynb` will need to include `Compare gene ranking` cells as seen in `human_analysis/` if you would like to compare genes/genes sets with some reference ranking. 
9. Update config file (see below)


The tables lists parameters required to run the analysis in this repository. These will need to be updated to run your own analysis. The * indicates optional parameters if you are comparing the ranks of your genes/gene sets with some reference ranking.

Note: Some of these parameters are required by the imported [ponyo](https://github.com/greenelab/ponyo) modules. 

| Name | Description |
| :--- | :---------- |
| local_dir| str: Parent directory on local machine to store intermediate results|
| dataset_name| str: Name for analysis directory, which contains the notebooks being run. For our analysis its named "human_analysis"|
| template_data_file | str: Path on your local machine where to write and store template gene expression data file|
| compendium_data_file | str: Path on your local machine where to write and store compendium gene expression data file|
| normalized_compendium_data_file | str: Path on your local machine where to write and store normalized compendium gene expression data file|
| shared_genes_file | str: Path on your local machine where to write and store genes that will be examined. These genes are the intersection of genes in your dataset versus a reference to ensure that there are not Nans in downstream analysis|
| scaler_transform_file | str: Path on your local machine where to write and store normalization transform to be used to process data for visualization|
| reference_gene_file* | str: Path to file that contains reference genes and their rank|
| refrence_gene_name_col| str: Name of the column header that contains the reference genes. This is found in reference_gene_file*|
| reference_rank_col | str: Name of the column header that contains the reference gene ranks. This is found in reference_gene_file*|
| NN_architecture | str: Name of neural network architecture to use. Format 'NN_<intermediate layer>_<latent layer>'|
| learning_rate| float: Step size used for gradient descent. In other words, it's how quickly the  methods is learning|
| batch_size | str: Training is performed in batches. So this determines the number of samples to consider at a given time|
| epochs | int: Number of times to train over the entire input dataset|
| kappa | float: How fast to linearly ramp up KL loss|
| intermediate_dim| int: Size of the hidden layer|
| latent_dim | int: Size of the bottleneck layer|
| epsilon_std | float: Standard deviation of Normal distribution to sample latent space|
| num_simulated| int: Simulate a compendia with these many experiments, created by shifting the template experiment these many times|
| project_id | str:  Experiment id to use as a template experiment|
| col_to_rank | str:  Name of column header from DE association statistic results. This column will be use to rank genes. Select `logFC`, `P.Value`, `adj.P.Val`, `t`|
| num_recount2_experiments | int:  Number of recount2 experiments to download. Note this will not be needed when we update the training to use all of recount2|
| compare_genes | bool:  1 if comparing gene ranks with reference gene ranks. 0 if just identifying generic genes and gene sets but not comparing against a reference.|