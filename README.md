# Detecting generic gene expression signals

**Alexandra J. Lee, Rani K. Powers, Dallas L. Mould, Dongbo Hu, Georgia Doing, Jake Crawford, James C. Costello, Deborah A. Hogan, Casey S. Greene**

**University of Pennsylvania, University of Colorado Anschutz Medical Campus, Dartmouth College**

**Rationale**: People performing differential expression (DE) analysis found that some genes and subsequent pathways are more likely to be differentially expressed even across a wide range of experimental designs ([Powers et. al., Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793 ); [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491)). 

Given that there exist these commonly DE genes and subsequent pathways, it is important to be able to distinguish between genes that are generic versus experiment or condition-specific. These specific genes may point to, say, those genes that are disease-specific and may reveal new insights into pathogenic mechanisms that might have gotten overlooked by examining all the DE genes in aggregate. And these disease-specific genes can also help to prioritize DEGs for follow-up wet experiments/functional experiments. For example, [Swindell et. al.](https://www.sciencedirect.com/science/article/pii/S0022202X16312465#fig3) identified IL-17A as an inducer of DEGs most uniquely elevated in psoriasis lesions compared to other skin diseases. Furthermore, clinical data demonstrating efficacy of anti-IL-17A therapy for moderate-to-severe psoriasis. In general being able to distinguish between generic vs context-specific signals is important to learning gene function and revealing insights into mechanism of disease.


**Challenge**: 
Current methods, including Powers et. al. and Crow et. al., to identify generic genes and pathways rely on manual curation. This curation effort included collecting a large set of samples with corresponding metadata, processing data and performing DE analysis to get ranked list of genes.

If you want to perform a new DE analysis in a different biological **context** (i.e. different organism, tissue, media) then you might not have the curated data available. Switching contexts will require re-curation. Similarly, using a different statistical method will require re-curation. This curation effort is very time intensive.


**Goal:**
To develop a method that can automatically distinguish between specific versus generic genes and pathways

**Results:**
We introduce a method based on latent space transformation in multi-layer neural networks that makes it possible to automate the analysis of generic genes, termed Specific cOntext Pattern Highlighting In Expression (SOPHIE). We validated that SOPHIE could recapitulate previously identified generic genes and pathways found using manual curation. These generic genes appear to act as gene hubs, which are associated with many biological processes. We applied SOPHIE to identify specific genes using a new experiment ---<TBD>

**Conclusions:**
We developed a method to automatically identify generic genes and pathways using public data without the need for curation. The generic signals identified from this method can be used to interpret study results and direct follow-up experiments.

## Directory Structure
| Folder/file | Description |
| --- | --- | 
| [configs](configs) | This folder contains configuration files used to set hyperparameters for the different experiments |
| [generic_expression_patterns_modules](generic_expression_patterns_modules) | This folder contains supporting functions that other notebooks in this repository will use |
| [human_cancer_analysis](human_cancer_analysis) | This folder contains analysis notebooks to validate generic signals using Powers et. al. dataset, which is composed of experiments testing the response of small molecule treatments in cancer cell lines, to train VAE |
| [human_general_analysis](human_general_analysis) | This folder contains analysis notebooks to validate generic signals using [recount2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6742427/) dataset, which contains a heterogeneous set of experiments, to train VAE |
| [multiplier_analysis](multiplier_analysis) | This folder contains analysis notebooks to coverage of generic genes across [MultiPLIER latent variables](https://www.cell.com/cell-systems/pdfExtended/S2405-4712\(19\)30119-X)) |
| [pseudomonas_analysis](pseudomonas_analysis) |  This folder contains analysis notebooks to identify specific and generic signals using *P. aeruginosa* dataset to train VAE |
| [new_experiment](new_experiment) |  This folder contains analysis notebooks to identify specific and generic signals using a new experiment and an existing VAE model|


## Usage

**How to run notebooks from generic-expression-patterns**

*Operating Systems:* Mac OS, Linux (Note: bioconda libraries not available in Windows)

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
6. Navigate to either the `pseudomonas_analysis`, `human_general_analysis` or `human_cancer_analysis` directories and run the notebooks in order.

*Note:* Running the `human_general_analysis/1_process_recount2_data.ipynb` notebook can take several days to run since the dataset is very large. If you would like to run only the analysis notebook (`human_general_analysis/2_identify_generic_genes_pathways.ipynb`) to generate the human analysis results found in the publication, you can update the config file to use the following file locations: 
* The normalized compendium data used for the analysis in the publication can be found [here](https://recount2.s3.amazonaws.com/normalized_recount2_compendium_data.tsv). 
* The Hallmark pathway database can be found [here](human_general_analysis/data/metadata/hallmark_DB.gmt)
* The processed template file can be found [here](human_general_analysis/data/processed_recount2_template.tsv)
* The scaler file can be found [here](human_general_analysis/data/scaler_transform_human.pickle)

**How to analyze your own data**

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
6.  Navigate to `new_experiment_example/find_specific_genes_in_new_experiment.ipynb` to see an example of how to run you analyze your own dataset using existing models
6. Create a configuration and metadata files for your analysis following the instructions in the `find_specific_genes_in_new_experiment.ipynb` notebook and the definitions below. Configuration files should be in `config/` directory. Metadata files should be within your analysis directory (`data/metadata/`). Here are the links to the compendium data needed:

<TO DO: Link to normalized and mapped compendium data on AWS>
7. Run notebook

*Note*:
* Your input dataset should be a matrix that is sample x gene
* The gene ids should be HGNC symbols (if using human data) or PA numbers (if using *P. aeruginosa* data)
* Your input dataset should be generated using the same platform as the model you plan to use (i.e. RNA-seq or array)
* Models available to use are: recount2 (human RNA-seq model found in `human_general_analysis/models`), Powers et. al. (human array model found in `human_cancer_analysis/models`), *P. aeruginosa* (*P. aeruginosa* array model found in `pseudomonas_analysis/models`)


The tables lists parameters required to run the analysis in this repository. These will need to be updated to run your own analysis. The * indicates optional parameters if you are comparing the ranks of your genes/gene sets with some reference ranking. The ** is only used if using `get_recount2_sra_subset` (in download_recount2_data.R).

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
| shared_genes_filename | str: Pickle file on your local machine where to write and store genes that will be examined. These genes are the intersection of genes in your dataset versus a reference to ensure that there are not Nans in downstream analysis|
| scaler_filename | str: Pickle file on your local machine where to write and store normalization transform to be used to process data for visualization|
| reference_gene_filename* | str: File that contains reference genes and their rank|
| reference_gene_name_col| str: Name of the column header that contains the reference genes. This is found in reference_gene_filename*|
| reference_rank_col | str: Name of the column header that contains the reference gene ranks. This is found in reference_gene_filename*|
| rank_genes_by | str:  Name of column header from DE association statistic results. This column will be use to rank genes. Select `logFC`, `P.Value`, `adj.P.Val`, `t` if using Limma. Select `log2FoldChange`, `pvalue`, `padj` if using DESeq.|
| pathway_DB_filename* | str: File that contains pathways to use for GSEA|
| gsea_statistic| str:  Statistic to use to rank genes for GSEA analysis. Select `logFC`, `P.Value`, `adj.P.Val`, `t` if using Limma. Select `log2FoldChange`, `pvalue`, `padj` if using DESeq.|
| rank_pathways_by | str:  Name of column header from GSEA association statistic results. This column will be use to rank pathways. Select `NES`, `padj` if using DESeq to rank genes.|
| NN_architecture | str: Name of neural network architecture to use. Format 'NN_<intermediate layer>_<latent layer>'|
| learning_rate| float: Step size used for gradient descent. In other words, it's how quickly the  methods is learning|
| batch_size | str: Training is performed in batches. So this determines the number of samples to consider at a given time|
| epochs | int: Number of times to train over the entire input dataset|
| kappa | float: How fast to linearly ramp up KL loss|
| intermediate_dim| int: Size of the hidden layer|
| latent_dim | int: Size of the bottleneck layer|
| epsilon_std | float: Standard deviation of Normal distribution to sample latent space|
| project_id | str:  Experiment id to use as a template experiment|
| count_threshold | int: Minimum count threshold to use to filter RNA-seq data|
| metadata_colname | str:  Header of experiment metadata file to indicate column containing sample ids. This is used to extract gene expression data associated with project_id|
| num_simulated| int: Simulate a compendia with these many experiments, created by shifting the template experiment these many times|
| num_recount2_experiments_to_download** | int:  Number of recount2 experiments to download. Note this will not be needed when we update the training to use all of recount2|

## Acknowledgements
We would like to thank David Nicholson, Ben Heil, Jake Crawford, Georgia Doing and Milton Pividori for insightful discussions and code review
