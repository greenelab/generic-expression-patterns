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

# # Plot
#
# This notebook plots the distribution of differences in mean ranks using 8 template experiments
#
# *Note: There were 10 template experiments but for one of the experiments DESeq didn't successfully run due to the distribution of data between samples being too similar

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import seaborn as sns
from ponyo import utils
import pickle

# ## Visualize the simulated data

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = "config_sophie_vs_trad.tsv"

params = utils.read_config(config_filename)

# +
# Load config params

# Local directory to store intermediate files
local_dir = params["local_dir"]

#
dataset_name = params["dataset_name"]

# File containing un-normalized template experiment
raw_template_filename = params["raw_template_filename"]

# ID for template experiment to be selected
project_id = params["project_id"]

# Identifier for simulated experiment
i = 0

# Pickle files saving specific and generic gene ids
template_specific_gene_ids_filename = params["template_specific_gene_ids_filename"]
generic_gene_ids_filename = "generic_gene_ids.pickle"
# -

# ## STOP

mean_diff_sophie = [
    18.279999999999973,
    157.5999999999999,
    -52.280000000000086,
    -73.39999999999998,
    -77.68999999999994,
    5.712121212121247,
    5.5,
    -121.57999999999993,
    187.05000000000007,
]

mean_diff_trad = [
    -146.63,
    97.46000000000004,
    -2.7000000000000455,
    -6.919999999999959,
    -92.82999999999993,
    -81.4353535353535,
    -192.5,
    -256.09000000000003,
    54.049999999999955,
]

# +
# Sample level
# -

"""mean_diff_sophie = [
    11.45999999999998,
    63.17999999999995,
    -38.05000000000007,
    -180.85000000000002,
    -20.860000000000014,
    52.52000000000004,
    47.93000000000001,
    -74.08000000000004,
]"""

"""mean_diff_trad = [
    -140.90999999999997,
    -34.6400000000001,
    -128.58000000000004,
    61.35000000000002,
    -125.88,
    -62.710000000000036,
    -262.26000000000005,
    -169.8,
]"""

"""mean_diff_sophie = [
    -11.45999999999998,
    -63.18000000000001,
    38.05000000000001,
    180.85000000000002,
    20.8599999999999,
    -16.639999999999986,
    -52.52000000000004,
    -47.930000000000064,
    74.08000000000004,
]"""

"""mean_diff_trad = [
    140.91000000000003,
    34.639999999999986,
    128.57999999999998,
    -61.349999999999966,
    125.88,
    -33.610000000000014,
    62.710000000000036,
    262.25999999999993,
    169.8,
]"""

# Make dataframe
# Each row in `mean_diff` is a different template experiment.
# The generic genes are the same across template experiment,
# but their rankings and therefore their mean may differ between template experiments.
mean_diff = pd.DataFrame(
    data={"sophie rank diff": mean_diff_sophie, "traditional rank diff": mean_diff_trad}
)

mean_diff

# Melt dataframe to use coloring in boxplot
mean_diff_melt = pd.melt(mean_diff)

"""# Plot coverage distribution given list of generic coverage, specific coverage
fig = sns.swarmplot(
    data=mean_diff_melt,
    x="variable",
    y="value",
    # notch=True,
    palette=["#2c7fb8", "lightgrey"],
    size=10,
)
fig.set_xlabel(None)
fig.set_xticklabels(["SOPHIE", "Traditional"], fontsize=14, fontname="Verdana")
fig.set_ylabel("Mean rank difference", fontsize=14, fontname="Verdana")
fig.tick_params(labelsize=14)
fig.set_title(
    "Mean rank difference between specific vs common DEGs",
    fontsize=16,
    fontname="Verdana",
)

fig.figure.savefig("validate_sophie_vs_trad.svg", format="svg", dpi=300)"""
