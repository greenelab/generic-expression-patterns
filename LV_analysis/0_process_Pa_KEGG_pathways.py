# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:generic_expression] *
#     language: python
#     name: conda-env-generic_expression-py
# ---

# # Process _P. aeruginosa_ KEGG pathways
#
# This notebook downloads _P. aeruginosa_ KEGG pathways from [here](https://raw.githubusercontent.com/greenelab/adage/master/Node_interpretation/pseudomonas_KEGG_terms.txt). These are the pathways used in the [ADAGE paper](https://msystems.asm.org/content/1/1/e00025-15) and formats it to be input to PLIER. PLIER expects pathway data to be a matrix of the form gene x pathway, where the values are 1 if the gene is contained with the pathway and 0 otherwise.
#
# This formated pathway data is used as in [PLIER R script](../generic_expression_patterns_modules/plier_util.R).

import pandas as pd

# Load data
expression_filename = "https://github.com/greenelab/adage/blob/2575a60804218db7f91402b955371bb60e5b00d6/Data_collection_processing/Pa_compendium_02.22.2014.pcl"
pa_kegg_pathway_filename = "https://github.com/greenelab/adage/blob/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

# Load expression data to get gene ids
expression_data = pd.read_csv(expression_filename, sep="\t", index_col=0, header=0)
expression_data.head()

# Load Pa KEGG pathway data
pa_kegg_pathway = pd.read_csv(
    pa_kegg_pathway_filename, sep="\t", index_col=0, header=None
)
pa_kegg_pathway.head()

# +
# Create new dataframe
# gene x pathway

gene_ids = list(expression_data.index)
pathway_names = list(pa_kegg_pathway.index)

out_pathway = pd.DataFrame(data=0, index=gene_ids, columns=pathway_names)

out_pathway

# +
# Fill in dataframe such that
# 1 if gene is contained within the pathway and 0 otherwise

for pathway_i in pathway_names:
    gene_ids_contained = pa_kegg_pathway.loc[pathway_i, 2].split(";")
    # Remove "." from gene ids
    gene_ids_contained = [gene_id.split(".")[0] for gene_id in gene_ids_contained]

    # Filter to gene ids that are contained within expression index
    gene_ids_contained = [
        gene_id for gene_id in gene_ids_contained if gene_id in gene_ids
    ]

    out_pathway.loc[gene_ids_contained, pathway_i] = 1

out_pathway.head()
# -

pathway_name_test = pa_kegg_pathway.index[0]
print(pathway_name_test)
assert out_pathway[pathway_name_test].sum() == len(
    pa_kegg_pathway.loc[pathway_name_test, 2].split(";")
), print(
    out_pathway[pathway_name_test].sum(),
    len(pa_kegg_pathway.loc[pathway_name_test, 2].split(";")),
)

# Manually check that counts match
out_pathway.sum().head(10)

pa_kegg_pathway.head(10)

# +
# Save dataframe
# This will be used in plier_util.R

out_pathway.to_csv("pa_kegg_pathway_processed.tsv", sep="\t")
