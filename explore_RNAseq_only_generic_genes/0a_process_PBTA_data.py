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

# # Process PBTA data

# +
# %load_ext autoreload
# %autoreload 2

import os
import pandas as pd
import pickle

from ponyo import utils
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# +
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
processed_template_filename = params["processed_template_filename"]
pbta_dir = os.path.join(local_dir, "openPBTA")
# -

# ## Load RDS objects

readRDS = ro.r["readRDS"]

polya_matrix = readRDS(
    os.path.join(pbta_dir, "pbta-gene-counts-rsem-expected_counts-collapsed.polya.rds")
)
ribo_matrix = readRDS(
    os.path.join(
        pbta_dir, "pbta-gene-counts-rsem-expected_counts-collapsed.stranded.rds"
    )
)

with localconverter(ro.default_converter + pandas2ri.converter):
    polya_matrix_values = ro.conversion.rpy2py(polya_matrix)
    ribo_matrix_values = ro.conversion.rpy2py(ribo_matrix)

polya_matrix_df = pd.DataFrame(
    data=polya_matrix_values,
    index=polya_matrix.rownames,
    columns=polya_matrix.colnames,
)
ribo_matrix_df = pd.DataFrame(
    data=ribo_matrix_values,
    index=ribo_matrix.rownames,
    columns=ribo_matrix.colnames,
)

print(polya_matrix_df.shape)
polya_matrix_df.head()

print(ribo_matrix_df.shape)
ribo_matrix_df.head()

# ## Get matching samples

# +
# Load metadata that maps RNA sample ids in expression matrices above
# to patient sample id
patient_metadata_filename = "https://raw.githubusercontent.com/kgaonkar6/OpenPBTA-analysis/532c29ab743bc643e687044bdb3e90241925186a/analyses/tp53_nf1_score/results/tp53_altered_status.tsv"

patient_metadata = pd.read_csv(
    patient_metadata_filename, sep="\t", index_col=0, header=0
)
# -

patient_metadata.head(10)

# +
# Select those patient sample ids (`sample_id) with multiple measurements
patient_metadata_tmp = patient_metadata[patient_metadata.index.value_counts() > 1]

# Select those with RNA sample ids (`Kids_First_Biospecimen_ID_RNA`) available
patient_metadata_selected = patient_metadata_tmp[
    patient_metadata_tmp["Kids_First_Biospecimen_ID_RNA"].isnull() == False
]

# +
# Create dictionary to map patient sample ids with those RNA ids
# for polyA-selection and ribo-depleted processing (column ids from gene expression matrices)
patient_sample_ids = list(patient_metadata_selected.index.unique())
polya_sample_ids = list(polya_matrix_df.columns)
ribo_sample_ids = list(ribo_matrix_df.columns)

patient_to_polya_id = {}
patient_to_ribo_id = {}
for patient_id in patient_sample_ids:
    rna_sample_ids = patient_metadata_selected.loc[
        patient_id, "Kids_First_Biospecimen_ID_RNA"
    ]
    for rna_sample_id in rna_sample_ids:
        if rna_sample_id in polya_sample_ids:
            patient_to_polya_id[patient_id] = rna_sample_id
        if rna_sample_id in ribo_sample_ids:
            patient_to_ribo_id[patient_id] = rna_sample_id

patient_to_polya_id
# -

patient_to_ribo_id

# +
# Select patient sample ids with both polyA-selected and ribo-depleted measurements
shared_patient_ids = list(patient_to_polya_id.keys() & patient_to_ribo_id.keys())

# Check that these patient ids were consistent with previous analysis comparing TP53 status across platform:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/930
shared_patient_ids
# -

# ## Select expression data

select_polya_ids = [patient_to_polya_id[x] for x in shared_patient_ids]
select_ribo_ids = [patient_to_ribo_id[x] for x in shared_patient_ids]

select_polya_expression = polya_matrix_df[select_polya_ids].T
select_ribo_expression = ribo_matrix_df[select_ribo_ids].T

# ## Format data matrix
#
# * Include only those genes that were used in our analysis
#     - Note: gene ENSEMBL ids already mapped to HGNC ids: https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq
# * Select only those samples with measurements from both polyA-selection and ribo-depleted protocols
# * Create metadata dataframe with grouping information for DE analysis.

# Read template file
template_SOPHIE = pd.read_csv(
    os.path.join(base_dir, "human_general_analysis", processed_template_filename),
    sep="\t",
    index_col=0,
    header=0,
)

template_SOPHIE.head()

# Get SOPHIE gene ids
SOPHIE_gene_ids = list(template_SOPHIE.columns)

# +
# Get shared gene ids between polyA and ribo
shared_platform_gene_ids = set(select_polya_expression.columns).intersection(
    select_ribo_expression.columns
)
print(len(shared_platform_gene_ids))

shared_gene_ids = list(set(shared_platform_gene_ids).intersection(SOPHIE_gene_ids))
print(len(shared_gene_ids))

# +
# Select shared genes
select_polya_expression = select_polya_expression[shared_gene_ids]
select_ribo_expression = select_ribo_expression[shared_gene_ids]

print(select_polya_expression.shape)
print(select_ribo_expression.shape)
# -

select_polya_expression.head()

select_ribo_expression.head()

# +
# Concatenate expression data
select_expression = pd.concat([select_polya_expression, select_ribo_expression])

select_expression.head(12)

# +
# Create metadata grouping matrix
polya_ids = list(select_polya_expression.index)
ribo_ids = list(select_ribo_expression.index)

sample_ids = polya_ids + ribo_ids
labels = [1] * 6 + [2] * 6

sample_grouping_metadata = pd.DataFrame(data={"Sample": sample_ids, "group": labels})

sample_grouping_metadata.set_index("Sample", inplace=True)
sample_grouping_metadata
# -

# ## Save

# +
expression_data_filename = "polya_ribo_expression.tsv"
sample_grouping_filename = "polya_ribo_sample_grouping.tsv"

select_expression.to_csv(expression_data_filename, sep="\t")
sample_grouping_metadata.to_csv(sample_grouping_filename, sep="\t")
