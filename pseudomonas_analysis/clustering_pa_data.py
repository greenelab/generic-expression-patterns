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

# # Clustering *P. aeruginosa* data
#
# This notebook will plot the clustering of P. aeruginosa data by strain type to demonstrate that clustering is perserved in VAE latent space.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import pandas as pd
import umap
import glob
import os
from keras.models import load_model
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

from numpy.random import seed

random_state = 1
seed(random_state)

# +
# Pa expression data
raw_pa_compendium_filename = "/home/alexandra/Documents/Data/Generic_expression_patterns/normalized_pseudomonas_compendium_data.tsv"

# Pa metadata
pa_metadata_filename = "data/metadata/sample_annotations.tsv"
# -

# Load data and metadata
raw_pa_compendium = pd.read_csv(
    raw_pa_compendium_filename, sep="\t", index_col=0, header=0
)
pa_metadata = pd.read_csv(pa_metadata_filename, sep="\t", index_col=0, header=0)

print(raw_pa_compendium.shape)
raw_pa_compendium.head()

print(pa_metadata.shape)
pa_metadata.head()

# ## Label expression samples

# Set "ml_data_source" (which corresponds to the sample ids in our expression matrix) as the index
pa_metadata.set_index("ml_data_source", inplace=True)

# Select and separate between samples that are using PAO1 strain and those using PA14 strain
pao1_sample_ids = pa_metadata.query("strain=='PAO1'").index.dropna()
pa14_sample_ids = pa_metadata.query("strain=='PA14'").index.dropna()

# Get shared sample ids
pao1_sample_ids_shared = set(raw_pa_compendium.index).intersection(pao1_sample_ids)
pa14_sample_ids_shared = set(raw_pa_compendium.index).intersection(pa14_sample_ids)

# Label samples
raw_pa_compendium["strain type"] = "other"
raw_pa_compendium.loc[pao1_sample_ids_shared, "strain type"] = "PAO1"
raw_pa_compendium.loc[pa14_sample_ids_shared, "strain type"] = "PA14"

raw_pa_compendium.head()

raw_pa_compendium_numeric = raw_pa_compendium.drop(columns=["strain type"])

# ## Plot

# +
# UMAP embedding of original input data

# pca = PCA(n_components=2)
# model = pca.fit(raw_pa_compendium_numeric)

# Get and save model
model = umap.UMAP(random_state=random_state).fit(raw_pa_compendium_numeric)

input_data_UMAPencoded = model.transform(raw_pa_compendium_numeric)
input_data_UMAPencoded_df = pd.DataFrame(
    data=input_data_UMAPencoded,
    index=raw_pa_compendium_numeric.index,
    columns=["1", "2"],
)
# Add label
input_data_UMAPencoded_df["strain type"] = raw_pa_compendium["strain type"]

input_data_UMAPencoded_df.head()
# -

# Colors
marker_colors = {
    "PA14": "#EF8B46",
    "PAO1": "#C6A9B5",
    "other": "#808080",
}

# +
# Plot
fig1 = ggplot(input_data_UMAPencoded_df, aes(x="1", y="2"))
fig1 += geom_point(aes(color="strain type"), alpha=0.2, size=3, stroke=0.8)
fig1 += scale_color_manual(values=marker_colors)
fig1 += labs(x="UMAP 1", y="UMAP 2", title="Pa compendium in gene space")
fig1 += theme_bw()
fig1 += theme(
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
fig1 += guides(colour=guide_legend(override_aes={"alpha": 1}))

print(fig1)

fig1.save("pa_raw_clustering.svg", format="svg", dpi=300)

# +
# Load model
vae_model_dir = "models/NN_2500_30/"

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
data_encoded = loaded_model.predict_on_batch(raw_pa_compendium_numeric)
data_encoded_df = pd.DataFrame(data_encoded, index=raw_pa_compendium_numeric.index)

# +
# UMAP embedding of VAE encoded data

# pca = PCA(n_components=2)
# model = pca.fit(data_encoded_df)
# Get and save model
model = umap.UMAP(random_state=random_state).fit(data_encoded_df)

input_data_UMAPencoded = model.transform(data_encoded_df)
input_data_UMAPencoded_df = pd.DataFrame(
    data=input_data_UMAPencoded, index=data_encoded_df.index, columns=["1", "2"]
)
# Add label
input_data_UMAPencoded_df["strain type"] = raw_pa_compendium["strain type"]

input_data_UMAPencoded_df.head()

# +
# Plot
fig2 = ggplot(input_data_UMAPencoded_df, aes(x="1", y="2"))
fig2 += geom_point(aes(color="strain type"), alpha=0.2, size=3, stroke=0.8)
fig2 += scale_color_manual(values=marker_colors)
fig2 += labs(x="UMAP 1", y="UMAP 2", title="Pa compendium in VAE space")
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

fig2.save("pa_vae_clustering.svg", format="svg", dpi=300)
