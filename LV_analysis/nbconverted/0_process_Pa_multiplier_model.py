# ---
# jupyter:
#   jupytext:
#     formats: LV_analysis//ipynb,LV_analysis/nbconverted//py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python [conda env:generic_expression] *
#     language: python
#     name: conda-env-generic_expression-py
# ---

# # Process _P. aeruginosa_ multiplier model

# +
# %load_ext autoreload
# %autoreload 2

import os
import pandas as pd
import pickle

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
# -

readRDS = ro.r["readRDS"]

# CHANGE LOCATION TO LOCAL WHEN COMMIT
multiplier_full_model = readRDS("Pa_compendium_PLIER_model.RDS")

# # Format multiplier Z data
#
# The Z data matrix contains the contribution (i.e. weight) per gene to each latent variable

multiplier_model_Z_matrix = multiplier_full_model.rx2("Z")

with localconverter(ro.default_converter + pandas2ri.converter):
    multiplier_model_Z_matrix_values = ro.conversion.rpy2py(multiplier_model_Z_matrix)

# +
column_header = [f"LV{i}" for i in range(1, 73)]

multiplier_model_Z_matrix_df = pd.DataFrame(
    data=multiplier_model_Z_matrix_values,
    index=multiplier_model_Z_matrix.rownames,
    columns=column_header,
)
# -

print(multiplier_model_Z_matrix_df.shape)
multiplier_model_Z_matrix_df.head()

# Save
multiplier_model_Z_matrix_df.to_csv("multiplier_Pa_model_z.tsv", sep="\t")

# # Format multiplier summary data
#
# This summary data matrix contains statistics about each LV - which pathways it was associated with and its significance score. This information is saved in the MultiPLIER model: https://github.com/greenelab/multi-plier/blob/7f4745847b45edf8fef3a49893843d9d40c258cf/23-explore_AAV_recount_LVs.Rmd

multiplier_model_matrix = multiplier_full_model.rx2("summary")

with localconverter(ro.default_converter + pandas2ri.converter):
    multiplier_model_matrix_values = ro.conversion.rpy2py(multiplier_model_matrix)

multiplier_model_matrix_df = pd.DataFrame(
    data=multiplier_model_matrix_values,
    index=multiplier_model_matrix.rownames,
    columns=multiplier_model_matrix.colnames,
)

multiplier_model_matrix_df.head()

# Save
multiplier_model_matrix_df.to_csv("multiplier_Pa_model_summary.tsv", sep="\t")
