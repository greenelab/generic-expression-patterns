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

# %load_ext autoreload
# %autoreload 2
import pandas as pd
import seaborn as sns

mean_diff_sophie = [
    -11.45999999999998,
    -63.18000000000001,
    38.05000000000001,
    180.85000000000002,
    20.8599999999999,
    -16.639999999999986,
    -52.52000000000004,
    -47.930000000000064,
    74.08000000000004,
]

mean_diff_trad = [
    140.91000000000003,
    34.639999999999986,
    128.57999999999998,
    -61.349999999999966,
    125.88,
    -33.610000000000014,
    62.710000000000036,
    262.25999999999993,
    169.8,
]

# Make dataframe
mean_diff = pd.DataFrame(
    data={"sophie rank diff": mean_diff_sophie, "traditional rank diff": mean_diff_trad}
)

mean_diff

# Melt dataframe to use coloring in boxplot
mean_diff_melt = pd.melt(mean_diff)

# Plot coverage distribution given list of generic coverage, specific coverage
fig = sns.swarmplot(
    data=mean_diff_melt,
    x="variable",
    y="value",
    # notch=True,
    palette=["#2c7fb8", "lightgrey"],
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
