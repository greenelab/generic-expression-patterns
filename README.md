# SOPHIE: Generative neural networks separate common and specific transcriptional responses

**Alexandra J. Lee, Dallas L. Mould, Jake Crawford, Dongbo Hu, Rani K. Powers, Georgia Doing, James C. Costello, Deborah A. Hogan, Casey S. Greene**

**University of Pennsylvania, University of Colorado Anschutz Medical Campus, Dartmouth College**

[![PDF Manuscript](https://img.shields.io/badge/manuscript-PDF-blue.svg)](https://www.biorxiv.org/content/10.1101/2021.05.24.445440v3)

There exist some genes and pathways that are differentially expressed across many gene expression experiments ([Powers et. al., Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793 ); [Crow et. al., PNAS 2019](https://www.pnas.org/content/116/13/6491)).
These common findings can obscure results that are specific to the context or experiment of interest, which are often what we hope to glean when using gene expression to generate mechanistic insights into cellular states and diseases.
Current methods, including Powers et. al. and Crow et. al., to identify common transcriptional signals rely on the manual curation and identical analysis of hundreds or thousands of additional experiments, which is inordinately time consuming and not a practical step in most analytical workflows.
If you want to perform a new DE analysis in a different biological **context** (i.e. different organism, tissue, media) then you might not have the curated data available. Switching contexts will require re-curation. Similarly, using a different statistical method will require re-curation.


We introduce a new approach to identify common patterns that uses generative neural networks to produce a null or background set of transcriptomic experiments.
Analyzing a target experiment against this automatically generated background set makes it straightforward to separate common and specific results.
This approach, called SOPHIE for Specific cOntext Pattern Highlighting In Expression data, can be applied to any new platform or species for which there is a large collection of unlabeled gene expression data.
Here, we apply SOPHIE to the analysis of both human and bacterial datasets, and use this method to highlight the ability to detect highly specific but low magnitude transcriptional signals that are biologically relevant.
The reusable notebooks for training neural networks and for the use of pre-trained generative models for the analysis of differential expression experiments may be broadly useful for the prioritization of specific findings in complex datasets.


## SOPHIE

This approach was named after one of the main characters from Hayao Miyazaki's animated film [Howl’s moving castle](https://en.wikipedia.org/wiki/Howl%27s_Moving_Castle_(film)).
Sophie’s outwardly appearance as an old woman despite being a young woman that has been cursed, demonstrates that the most obvious thing you see isn't always the truth.
This is the idea behind our approach, which allows users to identify specific gene expression signatures that can be masked by common background patterns.

SOPHIE applies [ponyo](https://github.com/greenelab/ponyo) to simulate gene expression experiments to use as a background set.
Then SOPHIE calculates differential expression statistics for each experiment in the background set.
The distribution of the log2 fold change of each gene across the simulated experiments can be used as a null to compare how changed a gene is in our target experiment of interest.
This approach allows investigators to distinguish common DEGs from context specific ones in their results.


## Directory Structure
| Folder | Description |
| --- | --- |
| [LV_analysis](LV_analysis) | This folder contains analysis notebooks to examine the potential role of common genes by looking at the coverage of common genes across [PLIER latent variables](https://www.cell.com/cell-systems/pdfExtended/S2405-4712\(19\)30119-X)) or [eADAGE latent variables](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5532071/), which are associated with known biological pathways.|
| [compare_experiments](compare_experiments) | This folder contains analysis notebooks to compare multiple SOPHIE results using the same template experiment and different template experiments. This analysis tests the robustness of SOPHIE results. |
| [configs](configs) | This folder contains configuration files used to set hyperparameters for the different experiments |
| [explore_RNAseq_only_generic_genes](explore_RNAseq_only_generic_genes) | This folder contains analysis notebooks testing different hypotheses to explain the subset of genes found to be common by SOPHIE trained on RNA-seq data but not found to be common in the manually curated array dataset. |
| [explore_data](explore_data) | This folder contains an analysis notebook visualizing the [recount2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6742427/) dataset to get a sense for the variation contained.|
| [expression_common_vs_other](expression_common_vs_other) | This folder contains an analysis notebook to determine if there are technical reasons that explain why common DEGs are common. Specifically this notebook is comparing the average expression of common genes versus other genes. |
| [figure_generation](figure_generation) | This folder contains a notebook to generate figures seen in the manuscript. |
| [generic_expression_patterns_modules](generic_expression_patterns_modules) | This folder contains supporting functions that other notebooks in this repository will use. |
| [human_cancer_analysis](human_cancer_analysis) | This folder contains analysis notebooks to validate common signals using [Powers et. al. dataset](https://academic.oup.com/bioinformatics/article/34/13/i555/5045793), which is composed of microarray experiments testing the response of small molecule treatments in human cancer cell lines, to train VAE. |
| [human_general_analysis](human_general_analysis) | This folder contains analysis notebooks to validate common signals using [recount2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6742427/) dataset, which contains a heterogeneous set of human RNA-seq experiments, to train VAE. |
| [human_general_array_analysis](human_general_array_analysis) | This folder contains analysis notebooks to validate common signals using [Crow et al](https://www.pnas.org/content/116/13/6491) dataset, which contains a heterogeneous set of human microarray experiments, to train VAE. |
| [network_analysis](network_analysis) |  This folder contains analysis notebooks to examine the potential role of common genes by looking at the clustering of generic genes within network communities.|
| [new_experiment](new_experiment) |  This folder contains analysis notebooks to identify specific and common signals using a new experiment and an existing VAE model|
| [new_model_experiment](new_model_experiment) |  This folder contains analysis notebooks to identify specific and common signals using a new experiment and training a new VAE model using your own compendium.|
| [other_enrichment_methods](other_enrichment_methods) |  This folder contains analysis notebooks to apply different gene set enrichment methods. The default method used is GSEA.|
| [pseudomonas_analysis](pseudomonas_analysis) |  This folder contains analysis notebooks to identify specific and common signals using *P. aeruginosa* microarray dataset to train VAE |
| [tests](tests) |  This folder contains notebooks to test the code in this repository. These notebooks run a small dataset across the analysis notebooks found in the `human_general_analysis` directory. |


## Usage

**How to reproduce the results and figures of the paper**

*Operating Systems:* Mac OS, Linux (Note: bioconda libraries not available in Windows)
*Note: *While the trends will be consistent, you may get slightly different resulting statistics and the plots may not look exactly the same as the paper.
Despite our best efforts to set seeds and version of python packages there is still some randomness we’re unable to control in the simulation process.

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
bash install.sh
```
This bash script uses the linux environment. If you are using a mac, you will need to update `install.sh` script to use `environment_mac.yml` and activate `generic_expression_mac`

6. Navigate to one of the analysis folders listed in the table above and run the jupyter notebooks in order.

*Note:* Running the `human_general_analysis/1_process_recount2_data.ipynb` notebook can take several days to run (this runtime was using a CPU) since the dataset is very large. If you would like to run only the analysis notebook (`human_general_analysis/2_identify_generic_genes_pathways.ipynb`) to generate the human analysis results found in the publication, you can update the config file to use the following file locations:
* The normalized compendium data used for the analysis in the publication can be found [here](https://storage.googleapis.com/recount2/normalized_recount2_compendium.tsv).
* The Hallmark pathway database can be found [here](human_general_analysis/data/metadata/hallmark_DB.gmt)
* The processed template file can be found [here](human_general_analysis/data/processed_recount2_template.tsv)
* The scaler file can be found [here](human_general_analysis/data/scaler_transform_human.pickle)

The runtime for training the VAE on the other datasets (`human_cancer_analysis/`, `human_general_array_analysis/`, `pseudomonas_analysis/`) were on the order of hours.
Here are links to the other stored large data files:

* normalized recount2 can be found [here](https://storage.googleapis.com/recount2/normalized_recount2_compendium.tsv)
* mapped recount2 can be found [here](https://storage.googleapis.com/recount2/mapped_recount2_compendium.tsv).
* normalized Powers et. al. can be found [here](https://storage.googleapis.com/powers_et_al/normalized_rani_compendium_filename.tsv).
* mapped Powers et. al. can be found [here](https://storage.googleapis.com/powers_et_al/mapped_rani_compendium.tsv).
* normalized _P. aeruginosa_ can be found [here](https://storage.googleapis.com/pseudomonas/normalized_pseudomonas_compendium_data.tsv).
* mapped _P. aeruginosa_ can be found [here](https://storage.googleapis.com/pseudomonas/processed_pseudomonas_compendium_data.tsv).

**How to analyze your own data using SOPHIE**

In order to run SOPHIE on your own analysis. Please visit the [sophie](https://github.com/greenelab/sophie) repository that includes all the scripts to run SOPHIE as well as templates documenting how to apply SOPHIE to your own dataset.

## Acknowledgements
We would like to thank David Nicholson, Ben Heil, Jake Crawford, Georgia Doing and Milton Pividori for insightful discussions and code review
