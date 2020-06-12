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

**Question**: Can we use our gene expression simulator to automatically generate null experiments for different contexts in order to overcome manual effort?

## Initial experiment
**Question**: Can our compendia simulation [link to github] identify the same common enriched pathways as filtered out in GSEA-InContext?

**Approach**:
1. Select treatment vs control experiment from recount2
2. Simulate 100 new experiments using experiment (1) as template
3. Perform DE analysis to get association statistics
4. Perform enrichment analysis using DEGs and rank pathways

**Hypothesis**: If we ranked pathways by the number of times that they showed up and we looked at those, we'd find the pathways that [GSEA-InContext](https://www.biorxiv.org/content/10.1101/659847v1) is designed to not find.

## Computational environment

Processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate generic_expression
```