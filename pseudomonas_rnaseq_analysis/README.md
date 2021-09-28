# Commonly differentially expressed genes

This analysis aims to using [SOPHIE](https://www.biorxiv.org/content/10.1101/2021.05.24.445440v1) to identify genes that are commonly differentially expressed in our PAO1 and PA14 compendia.
These commonly differentially expressed genes will be compared to the core and accessory genes to examine the different transcriptional behvior of core versus accessory genes.

## Dependencies

We use Python and Jupyter notebooks to load, process, and analyze the data.
Since the analyses are using the SOPHIE framework, we manually copied the scripts (`generic_expression_patterns_modules`) and environment files from the [generic-expression-patterns](https://github.com/greenelab/generic-expression-patterns) repository.

The environment dependencies are specified in the `environment.yml` file in this directory, and can be installed by running the following terminal commands from this directory:

```bash
conda env create --file environment.yml
conda activate generic_expression
```
