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

# # Visualize gene expression
#
# This notebook visualizes the gene expression data for the template and simulated experiments in order to:
# 1. Validate that the structure of the gene expression data and simulated data are consistent
# 2. To visualize the signal that is in the experiments

# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

# +
import os
import pandas as pd
import umap
import pickle
import glob
import seaborn as sns
from sklearn.decomposition import PCA
from keras.models import load_model
import plotnine as pn

from ponyo import utils
from generic_expression_patterns_modules import plot
# -

# ## Load config parameters

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)

params = utils.read_config(config_filename)

# +
# Load config params

local_dir = params["local_dir"]
project_id = params["project_id"]
num_simulated = params["num_simulated"]

pval_name = "adj.P.Val"
logFC_name = "logFC"
run = 0

# Manual settings to visualize/troubleshoot volcano plots for other datasets
# Will pull these out to archive later

"""vae_model_dir = params['vae_model_dir']
template_filename = params['mapped_template_filename']
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = params['scaler_filename']"""


# Settings for running visualization using pseudomonas config file
vae_model_dir = os.path.join(base_dir, "pseudomonas_analysis", "models", "NN_2500_30")
template_filename = params["processed_template_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]
scaler_filename = params["scaler_filename"]

"""# Settings for running visualization using human cancer config file
vae_model_dir = os.path.join(base_dir,"human_cancer_analysis", "models", "NN_2500_30")
template_filename = params['processed_template_filename']
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = params['scaler_filename']"""


"""# Settings for running visualization using human_general config file
vae_model_dir = os.path.join(base_dir,"human_general_analysis", "models", "NN_2500_30")
template_filename = os.path.join(base_dir,"human_general_analysis", params['processed_template_filename'])
normalized_compendium_filename = params['normalized_compendium_filename']
scaler_filename = os.path.join(base_dir, "human_general_analysis", params['scaler_filename'])"""
# -

# ## Volcano plots

# +
# Check number of DEGs
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)

template_DE_stats = pd.read_csv(
    template_DE_stats_filename, sep="\t", header=0, index_col=0
)

selected = template_DE_stats[
    (template_DE_stats[pval_name] < 0.01) & (abs(template_DE_stats[logFC_name]) > 1)
]
print(selected.shape)
# -

plot.make_volcano_plot_template(
    template_DE_stats_filename, project_id, pval_name, logFC_name
)

# +
simulated_DE_stats_dir = os.path.join(local_dir, "DE_stats")

plot.make_volcano_plot_simulated(
    simulated_DE_stats_dir,
    project_id,
    pval_name,
    logFC_name,
    num_simulated,
    5,
    5,
    20,
    15,
)
# -

# ## Plot distribution of DE stats

sns.distplot(template_DE_stats[logFC_name], kde=False)

# +
simulated_DE_stats_filename = os.path.join(
    simulated_DE_stats_dir, f"DE_stats_simulated_data_{project_id}_{run}.txt"
)

simulated_DE_stats = pd.read_csv(
    simulated_DE_stats_filename, sep="\t", header=0, index_col=0
)

sns.distplot(simulated_DE_stats[logFC_name], kde=False)
# -

# ## PCA of latent space

# Get decoded simulated experiment
simulated_filename = os.path.join(
    local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_{run}.txt"
)

normalized_compendium_data = pd.read_csv(
    normalized_compendium_filename, sep="\t", index_col=0, header=0
)
template_data = pd.read_csv(template_filename, sep="\t", index_col=0, header=0)
simulated_data = pd.read_csv(simulated_filename, sep="\t", index_col=0, header=0)

print(template_data.shape)
template_data

sns.distplot(template_data.iloc[0:2].mean(), kde=False)
sns.distplot(template_data.iloc[2:].mean(), kde=False)

sns.distplot(simulated_data.iloc[0:2].mean(), kde=False)
sns.distplot(simulated_data.iloc[2:].mean(), kde=False)

# +
# Normalize simulated_data
# Load pickled file
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)

normalized_simulated_data = scaler.transform(simulated_data)

normalized_simulated_data = pd.DataFrame(
    normalized_simulated_data,
    columns=simulated_data.columns,
    index=simulated_data.index,
)

print(normalized_simulated_data.shape)
normalized_simulated_data.head()

# +
# If template experiment included in training compendium
# Get normalized template data
sample_ids = list(template_data.index)
normalized_template_data = normalized_compendium_data.loc[sample_ids]

print(normalized_template_data.shape)
normalized_template_data.head()
# -

"""# If template experiment NOT included in training compendium
with open(scaler_filename, "rb") as scaler_fh:
    scaler = pickle.load(scaler_fh)

normalized_template_data = scaler.transform(template_data)

normalized_template_data = pd.DataFrame(
    normalized_template_data,
    columns=template_data.columns,
    index=template_data.index,
)"""

# Label samples
normalized_compendium_data["sample group"] = "compendium"
normalized_template_data["sample group"] = "template"
normalized_simulated_data["sample group"] = "simulated"

normalized_all_data = pd.concat(
    [normalized_template_data, normalized_simulated_data, normalized_compendium_data]
)

# +
# Plot

# Drop label column
normalized_all_data_numeric = normalized_all_data.drop(["sample group"], axis=1)

model = umap.UMAP(random_state=1).fit(normalized_all_data_numeric)

normalized_all_data_UMAPencoded = model.transform(normalized_all_data_numeric)
normalized_all_data_UMAPencoded_df = pd.DataFrame(
    data=normalized_all_data_UMAPencoded,
    index=normalized_all_data.index,
    columns=["1", "2"],
)

# Add back label column
normalized_all_data_UMAPencoded_df["sample group"] = normalized_all_data["sample group"]

# Plot
fig = pn.ggplot(normalized_all_data_UMAPencoded_df, pn.aes(x="1", y="2"))
fig += pn.geom_point(pn.aes(color="sample group"), alpha=0.4)
fig += pn.labs(x="UMAP 1", y="UMAP 2", title="Gene expression data in gene space")
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)
fig += pn.scale_color_manual(["#bdbdbd", "red", "blue"])
fig += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))

print(fig)
# -

# ## PCA in latent space

# Model files
model_encoder_filename = glob.glob(os.path.join(vae_model_dir, "*_encoder_model.h5"))[0]
weights_encoder_filename = glob.glob(
    os.path.join(vae_model_dir, "*_encoder_weights.h5")
)[0]
model_decoder_filename = glob.glob(os.path.join(vae_model_dir, "*_decoder_model.h5"))[0]
weights_decoder_filename = glob.glob(
    os.path.join(vae_model_dir, "*_decoder_weights.h5")
)[0]

# +
# Load saved models
loaded_model = load_model(model_encoder_filename, compile=False)
loaded_decode_model = load_model(model_decoder_filename, compile=False)

loaded_model.load_weights(weights_encoder_filename)
loaded_decode_model.load_weights(weights_decoder_filename)
# -

# PCA model
pca = PCA(n_components=2)

# +
# Encode compendium
normalized_compendium = normalized_compendium_data.drop(["sample group"], axis=1)
compendium_encoded = loaded_model.predict_on_batch(normalized_compendium)

compendium_encoded_df = pd.DataFrame(
    data=compendium_encoded, index=normalized_compendium.index
)

# Get and save PCA model
model1 = pca.fit(compendium_encoded_df)

compendium_PCAencoded = model1.transform(compendium_encoded_df)

compendium_PCAencoded_df = pd.DataFrame(
    data=compendium_PCAencoded, index=compendium_encoded_df.index, columns=["1", "2"]
)

# Add label
compendium_PCAencoded_df["sample group"] = "compendium"

# +
# Encode template experiment
normalized_template_data = normalized_template_data.drop(["sample group"], axis=1)

template_encoded = loaded_model.predict_on_batch(normalized_template_data)
template_encoded_df = pd.DataFrame(
    data=template_encoded, index=normalized_template_data.index
)

template_PCAencoded = model1.transform(template_encoded_df)

template_PCAencoded_df = pd.DataFrame(
    data=template_PCAencoded, index=template_encoded_df.index, columns=["1", "2"]
)

# Add back label column
template_PCAencoded_df["sample group"] = "template"

# +
# Use stored encoded simulated data
# Note: We cannot encode the decoded simulated experiment since we are not using tied weights
# Re-encoded the decoded simulated experiment will not yield a linear latent space shift
encoded_simulated_filename = os.path.join(
    local_dir,
    "pseudo_experiment",
    f"selected_simulated_encoded_data_{project_id}_{run}.txt",
)

simulated_encoded_df = pd.read_csv(
    encoded_simulated_filename, header=0, sep="\t", index_col=0
)

sample_ids = list(template_data.index)
simulated_encoded_df = simulated_encoded_df.loc[sample_ids]

simulated_PCAencoded = model1.transform(simulated_encoded_df)

simulated_PCAencoded_df = pd.DataFrame(
    data=simulated_PCAencoded, index=simulated_encoded_df.index, columns=["1", "2"]
)

# Add back label column
simulated_PCAencoded_df["sample group"] = "simulated"

# +
# Concatenate dataframes
combined_PCAencoded_df = pd.concat(
    [compendium_PCAencoded_df, template_PCAencoded_df, simulated_PCAencoded_df]
)

print(combined_PCAencoded_df.shape)
combined_PCAencoded_df.head()

# +
# Plot
fig1 = pn.ggplot(combined_PCAencoded_df, pn.aes(x="1", y="2"))
fig1 += pn.geom_point(pn.aes(color="sample group"), alpha=0.4)
fig1 += pn.labs(x="PC 1", y="PC 2", title="Gene expression data in latent space")
fig1 += pn.theme_bw()
fig1 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)
fig1 += pn.scale_color_manual(["#bdbdbd", "red", "blue"])
fig1 += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))

print(fig1)
# -

# ## UMAP of latent space

# +
# Get and save PCA model
model2 = umap.UMAP(random_state=1).fit(compendium_encoded_df)

compendium_UMAPencoded = model2.transform(compendium_encoded_df)

compendium_UMAPencoded_df = pd.DataFrame(
    data=compendium_UMAPencoded, index=compendium_encoded_df.index, columns=["1", "2"]
)

# Add label
compendium_UMAPencoded_df["sample group"] = "compendium"

# +
template_UMAPencoded = model2.transform(template_encoded_df)

template_UMAPencoded_df = pd.DataFrame(
    data=template_UMAPencoded, index=template_encoded_df.index, columns=["1", "2"]
)

# Add back label column
template_UMAPencoded_df["sample group"] = "template"

# +
simulated_UMAPencoded = model2.transform(simulated_encoded_df)

simulated_UMAPencoded_df = pd.DataFrame(
    data=simulated_UMAPencoded, index=simulated_encoded_df.index, columns=["1", "2"]
)

# Add back label column
simulated_UMAPencoded_df["sample group"] = "simulated"

# +
# Concatenate dataframes
combined_UMAPencoded_df = pd.concat(
    [compendium_UMAPencoded_df, template_UMAPencoded_df, simulated_UMAPencoded_df]
)

print(combined_UMAPencoded_df.shape)
combined_UMAPencoded_df.head()

# +
# Plot
fig2 = pn.ggplot(combined_UMAPencoded_df, pn.aes(x="1", y="2"))
fig2 += pn.geom_point(pn.aes(color="sample group"), alpha=0.4)
fig2 += pn.labs(x="UMAP 1", y="UMAP 2", title="Gene expression data in latent space")
fig2 += pn.theme_bw()
fig2 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)
fig2 += pn.scale_color_manual(["#bdbdbd", "red", "blue"])
fig2 += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))

print(fig2)
