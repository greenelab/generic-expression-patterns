# Generic genes network analysis

This analysis aims to explore characteristics of generic genes in the context of network representations of *P. aeruginosa* gene similarity/connectedness.

Specifically, in these analyses, we are using a network generated using correlation values between genes, where these correlations are defined based on eADAGE weight vectors. This network construction process is described in more detail in the [ADAGE signature analysis manuscript](https://doi.org/10.1186/s12859-017-1905-4) (see sections on "Signature interpretation" and "GeneNetwork" analysis).

## Dependencies

We use Python and Jupyter notebooks to load, process, and visualize network data. Since these analyses have several network analysis-specific dependencies, including the [graph-tool](https://graph-tool.skewed.de/) and  [python-igraph](https://igraph.org/python/) packages, they use a different Conda environment from the rest of the `generic-expression-patterns` repo.

The environment dependencies are specified in the `environment.yml` file in this directory, and can be installed by running the following terminal commands from this directory:

```bash
conda env create --file environment.yml
conda activate generic_network
```
