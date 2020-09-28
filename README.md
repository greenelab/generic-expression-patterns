# Rank pathways analysis

**Alexandra J Lee, James C Costello and Casey S Greene**

**May 2020**

**University of Pennsylvania, University of Colorado Anschutz Medical Campus**

**Rationale**: People performing differential expression (DE) analysis find that some genes and subsequent pathways are more likely to be differentially expressed even across a wide range of experimental designs ([Powers et. al., Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793 ); [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491)). 

Given that there exists these commonly DE genes and subsequent pathways, it is important to be able to distinguish between genes that are generic versus experiment or condition-specific. These specific genes may point to, say, those genes that are disease-specific and may reveal new insights into pathogenic mechanisms that might have gotten overlooked by examining all the DE genes in aggregate. And these disease-specific genes can also help to prioritize DEGs for follow-up wet experiments/functional experiments. For example, [Swindell et. al.](https://www.sciencedirect.com/science/article/pii/S0022202X16312465#fig3) identified IL-17A as an inducer of DEGs most uniquely elevated in psoriasis lesions compared to other skin diseases. Furthermore, clinical data demonstrating efficacy of anti-IL-17A therapy for moderate-to-severe psoriasis. In general being able to distinguish between generic vs context-specific signals is important to learning gene function and revealing insights into mechanism of disease.


**Challenge**: 
Current methods, including Powers et. al. and Crow et. al., to identify generic genes and pathways rely on manual curation. This curation effort included collecting a large set of samples with corresponding metadata, process data and perform DE analysis to get ranked list of genes.

If you want to perform a new DE analysis in a different biological **context** (i.e. different organism, tissue, media) then you might not have the curated data available. Switching contexts will require a lot of manual effort. Similarly, using a different statistical method will require re-curation effort


**Goal of this study:**
* To show that our compendia simulation method, [ponyo](https://github.com/greenelab/ponyo) can automatically identify generic genes and pathways

**Results:**
* We found a set of general generic genes (i.e. genes found to be generic in both recount2 and crow et. al., which contain a mix of experiments)
* We developed a method to automatically identify generic genes in different contexts without having to perform experiments and curate. 


## How to run notebooks from generic-expression-patterns

**Operating Systems:** Mac OS, Linux

In order to run this simulation on your own gene expression data the following steps should be performed:

First you need to set up your local repository: 
1. Download and install [github's large file tracker](https://git-lfs.github.com/).
2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
3. Clone the `generic-expression-patterns` repository by running the following command in the terminal:
```
git clone https://github.com/greenelab/generic-expression-patterns.git
```
Note: Git automatically detects the LFS-tracked files and clones them via http.
4. Navigate into cloned repo by running the following command in the terminal:
```
cd generic-expression-patterns
```
5. Set up conda environment by running the following command in the terminal:
```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate generic_expression

pip install -e .
```
6. Navigate to either the `pseudomonas_analysis` or `human_analysis` directories and run the notebooks in order.

*Note:* Running the `human_analysis/1_process_recount2_data.ipynb` notebook can take a while since the dataset is very large. If you would like to run only the analysis  (`2_identify_generic_genes_pathways.ipynb`) to generate the human analysis results found in the publication, you can update the config file to use the following file locations: 
* The normalized compendium data used for the analysis in the publication can be found [here](https://recount2.s3.amazonaws.com/normalized_recount2_compendium_data.tsv). 
* The mapped template file can be found ____
* The scaler file can be found ____

## How to run using your own data

In order to run this simulation on your own gene expression data the following steps should be performed:

First you need to set up your local repository: 
1. Download and install [github's large file tracker](https://git-lfs.github.com/).
2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
3. Clone the `generic-expression-patterns` repository by running the following command in the terminal:
```
git clone https://github.com/greenelab/generic-expression-patterns.git
```
Note: Git automatically detects the LFS-tracked files and clones them via http.
4. Navigate into cloned repo by running the following command in the terminal:
```
cd generic-expression-patterns
```
5. Set up conda environment by running the following command in the terminal:
```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate generic_expression

pip install -e .
```
5. Create a new analysis folder in the main directory. This is equivalent to the `human_analysis` directory
6. Copy jupyter notebooks (1_process_data.ipynb, 2_identify_generic_genes_pathways.ipynb) into analysis directory.
7. Customize `1_process_data.ipynb` to generate the following saved files: 1) compendium of gene expression data, 2) template experiment data, 3) normalized gene expression compendium. See examples of this processing in `human_analysis/` and `pseudomonas_analysis/`.
8. Required data files:
* Gene expression compendium to use as input data. Specify path of data file in config file.
* Metadata matrix (sample x experimental metadata) with sample id as index. Add file to `analysis_dir/data/metadata/`. Specify path of metadata file in config file.
* Sample grouping matrix (sample x group) with group = [1 if control; 2 if case]. Add file to `analysis_dir/data/metadata/<project_id>_groups.tsv`.
9. `2_identify_generic_genes_pathways.ipynb` will need to include `Compare gene ranking` cells as seen in `human_analysis/` if you would like to compare genes/genes sets with some reference ranking. 
10. Update config file (see below)


The tables lists parameters required to run the analysis in this repository. These will need to be updated to run your own analysis. The * indicates optional parameters if you are comparing the ranks of your genes/gene sets with some reference ranking.

Note: Some of these parameters are required by the imported [ponyo](https://github.com/greenelab/ponyo) modules. 

| Name | Description |
| :--- | :---------- |
| local_dir| str: Parent directory on local machine to store intermediate results|
| dataset_name| str: Name for analysis directory, which contains the notebooks being run. For our analysis its named "human_analysis"|
| raw_template_filename | str: Downloaded template gene expression data file|
| mapped_template_filename | str: Template gene expression data file after replacing gene ids in header|
| processed_template_filename | str: Template gene expression data file after removing samples and genes|
| raw_compendium_filename | str: Downloaded compendium gene expression data file|
| mapped_compendium_filename | str: Compendium gene expression data file after replacing gene ids in header|
| normalized_compendium_filename | str: Normalized compendium gene expression data file|
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
| rank_genes_by | str:  Name of column header from DE association statistic results. This column will be use to rank genes. Select `logFC`, `P.Value`, `adj.P.Val`, `t` if using Limma. Select `log2FoldChange`, `pvalue`, `padj` if using DESeq.|
| rank_pathways_by | str:  Name of column header from GSEA association statistic results. This column will be use to rank pathways. Select `NES`, `padj` if using DESeq to rank genes.|
| num_recount2_experiments_to_download | int:  Number of recount2 experiments to download. Note this will not be needed when we update the training to use all of recount2|
| gsea_statistic| str:  Statistic to use to rank genes for GSEA analysis. Select `logFC`, `P.Value`, `adj.P.Val`, `t` if using Limma. Select `log2FoldChange`, `pvalue`, `padj` if using DESeq.|
| compare_genes | bool:  1 if comparing gene ranks with reference gene ranks. 0 if just identifying generic genes and gene sets but not comparing against a reference.|
