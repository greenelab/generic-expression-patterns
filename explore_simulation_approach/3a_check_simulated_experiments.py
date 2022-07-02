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

# +
# %load_ext autoreload
# %autoreload 2
# %load_ext rpy2.ipython
# %matplotlib inline
import os
import glob
import pickle
import umap
import pandas as pd
import seaborn as sns
from keras.models import load_model
from ponyo import utils
from sklearn.decomposition import PCA

from plotnine import (
    ggplot,
    labs,
    geom_line,
    geom_point,
    geom_errorbar,
    aes,
    ggsave,
    theme_bw,
    theme,
    xlim,
    ylim,
    facet_wrap,
    scale_color_manual,
    guides,
    guide_legend,
    element_blank,
    element_text,
    element_rect,
    element_line,
    coords,
    options,
)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = "config_sophie_vs_trad.tsv"

params = utils.read_config(config_filename)

# +
# Load config params

# Local directory to store intermediate files
local_dir = params["local_dir"]

# ID for template experiment to be selected
project_id = params["project_id"]

# File containing un-normalized template experiment
raw_template_filename = params["raw_template_filename"]

# Normalized compendium filename
normalized_compendium_filename = params["normalized_compendium_filename"]

# Directory containing trained VAE model
vae_model_dir = params["vae_model_dir"]

# Pickle files saving specific and generic gene ids
template_specific_gene_ids_filename = params["template_specific_gene_ids_filename"]
generic_gene_ids_filename = "generic_gene_ids.pickle"

# Identifier for simulated experiment
i = 7

# File to scaler transform data to 0-1 range
scaler_transform_filename = params["scaler_filename"]
# -

# ## Check clustering

simulated_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_data_{project_id}_{i}.txt",
)

simulated_experiment = pd.read_csv(simulated_filename, sep="\t", index_col=0, header=0)

print(simulated_experiment.shape)
simulated_experiment.head()

# +
# Get common genes
# Get specific genes

# Load pickled file
with open(template_specific_gene_ids_filename, "rb") as specific_fh:
    specific_gene_ids = pickle.load(specific_fh)

with open(generic_gene_ids_filename, "rb") as generic_fh:
    generic_gene_ids = pickle.load(generic_fh)
# -

# Get NA genes
all_gene_ids = simulated_experiment.columns
all_gene_ids_tmp = all_gene_ids.difference(specific_gene_ids)
na_gene_ids = all_gene_ids_tmp.difference(generic_gene_ids)

# Simulated data subsets
simulated_specific_df = simulated_experiment[specific_gene_ids]
simulated_common_df = simulated_experiment[generic_gene_ids]
simulated_na_df = simulated_experiment[na_gene_ids]

print(simulated_specific_df.shape)
simulated_specific_df

print(simulated_common_df.shape)
simulated_common_df

print(simulated_na_df.shape)
simulated_na_df

f = sns.clustermap(simulated_specific_df.T, cmap="viridis")
f.fig.suptitle("Simulated experiment specific genes")

f = sns.clustermap(simulated_common_df.T, cmap="viridis")
f.fig.suptitle("Simulated experiment common genes")

f = sns.clustermap(simulated_na_df.T, cmap="viridis")
f.fig.suptitle("Simulated experiment NA genes")

# ## Check latent space distribution

normalized_compendium = pd.read_csv(
    normalized_compendium_filename, sep="\t", index_col=0, header=0
)

# +
# Files
model_encoder_filename = glob.glob(os.path.join(vae_model_dir, "*_encoder_model.h5"))[0]

weights_encoder_filename = glob.glob(
    os.path.join(vae_model_dir, "*_encoder_weights.h5")
)[0]

model_decoder_filename = glob.glob(os.path.join(vae_model_dir, "*_decoder_model.h5"))[0]

weights_decoder_filename = glob.glob(
    os.path.join(vae_model_dir, "*_decoder_weights.h5")
)[0]

# Load saved models
loaded_model = load_model(model_encoder_filename)
loaded_decode_model = load_model(model_decoder_filename)

loaded_model.load_weights(weights_encoder_filename)
loaded_decode_model.load_weights(weights_decoder_filename)
# -

# Encode normalized expression data
compendium_encoded = loaded_model.predict_on_batch(normalized_compendium)
compendium_encoded_df = pd.DataFrame(
    compendium_encoded, index=normalized_compendium.index
)

# +
# PCA embedding of VAE encoded data

random_state = 1
pca = PCA(n_components=2)
model = pca.fit(compendium_encoded_df)

compendium_data_PCencoded = model.transform(compendium_encoded_df)
compendium_data_PCencoded_df = pd.DataFrame(
    data=compendium_data_PCencoded,
    index=compendium_encoded_df.index,
    columns=["1", "2"],
)

print(compendium_data_PCencoded_df.shape)
compendium_data_PCencoded_df.head()

# +
# Plot
fig1 = ggplot(compendium_data_PCencoded_df, aes(x="1", y="2"))
fig1 += geom_point(alpha=0.2, size=3, stroke=0.8)
fig1 += labs(x="PC 1", y="PC 2", title="Background compendium in VAE space")
fig1 += theme_bw()
fig1 += theme(
    legend_title_align="center",
    plot_background=element_rect(fill="white"),
    legend_key=element_rect(fill="white", colour="white"),
    legend_title=element_text(family="sans-serif", size=15),
    legend_text=element_text(family="sans-serif", size=12),
    plot_title=element_text(family="sans-serif", size=15),
    axis_text=element_text(family="sans-serif", size=12),
    axis_title=element_text(family="sans-serif", size=15),
)
fig1 += guides(colour=guide_legend(override_aes={"alpha": 1}))

print(fig1)
# -

# What does the template compendium look like in the encoded latent space? Perhaps the encoder is compressing the difference in the perturbed and control samples too much.
#
# If that is the case, then trying to increase the variance between the compendium experiments might help reduce the compression in the latent space.

# +
# File storing normalized template experiment
normalized_template_filename = params["normalized_template_filename"]

normalized_template = pd.read_csv(
    normalized_template_filename, sep="\t", index_col=0, header=0
)
# -

normalized_template.head()

# Encode template experiment into VAE space
template_data_encoded = loaded_model.predict_on_batch(normalized_template)
template_data_encoded_df = pd.DataFrame(
    template_data_encoded, index=normalized_template.index
)

# +
simulated_encoded_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_encoded_data_{project_id}_{i}.txt",
)

simulated_data_encoded_df = pd.read_csv(
    simulated_encoded_filename, sep="\t", index_col=0, header=0
)

# +
# Encode into PC space
template_data_PCencoded = model.transform(template_data_encoded_df)
template_data_PCencoded_df = pd.DataFrame(
    data=template_data_PCencoded,
    index=template_data_encoded_df.index,
    columns=["1", "2"],
)

print(template_data_PCencoded_df.shape)
template_data_PCencoded_df.head(10)

# +
# Encode into PC space
simulated_data_PCencoded = model.transform(simulated_data_encoded_df)
simulated_data_PCencoded_df = pd.DataFrame(
    data=simulated_data_PCencoded,
    index=simulated_data_encoded_df.index,
    columns=["1", "2"],
)

print(simulated_data_PCencoded_df.shape)
simulated_data_PCencoded_df.head(10)

# +
# Add labels
compendium_data_PCencoded_df["label"] = "compendium background"
template_data_PCencoded_df.loc[
    normalized_template.index.str.contains("perturb"), "label"
] = "template perturb samples"
template_data_PCencoded_df.loc[
    normalized_template.index.str.contains("control"), "label"
] = "template control samples"

simulated_data_PCencoded_df.loc[
    simulated_experiment.index.str.contains("perturb"), "label"
] = "simulated perturb samples"
simulated_data_PCencoded_df.loc[
    simulated_experiment.index.str.contains("control"), "label"
] = "simulated control samples"

# +
all_data_PCencoded_df = pd.concat(
    [
        compendium_data_PCencoded_df,
        template_data_PCencoded_df,
        simulated_data_PCencoded_df,
    ]
)

print(all_data_PCencoded_df.shape)
all_data_PCencoded_df.head(20)
# -

# Colors
marker_colors = {
    "template perturb samples": "red",
    "template control samples": "blue",
    "simulated perturb samples": "magenta",
    "simulated control samples": "cyan",
    "compendium background": "grey",
}

# +
# Plot
fig2 = ggplot(all_data_PCencoded_df, aes(x="1", y="2"))
fig2 += geom_point(aes(color="label"), alpha=0.1, size=3, stroke=0.8)
fig2 += scale_color_manual(values=marker_colors)
fig2 += labs(x="PC 1", y="PC 2", title="Compendium in VAE space")
fig2 += theme_bw()
fig2 += theme(
    # figure_size=(6, 10),
    legend_title_align="center",
    plot_background=element_rect(fill="white"),
    legend_key=element_rect(fill="white", colour="white"),
    legend_title=element_text(family="sans-serif", size=15),
    legend_text=element_text(family="sans-serif", size=12),
    plot_title=element_text(family="sans-serif", size=15),
    axis_text=element_text(family="sans-serif", size=12),
    axis_title=element_text(family="sans-serif", size=15),
)
fig2 += guides(colour=guide_legend(override_aes={"alpha": 1}))

print(fig2)
