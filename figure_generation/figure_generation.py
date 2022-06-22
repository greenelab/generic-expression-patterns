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
#     display_name: Python [conda env:generic_expression_new]
#     language: python
#     name: conda-env-generic_expression_new-py
# ---

# # Figure generation

from IPython.display import Image, display, SVG
import svgutils.transform as sg
import numpy as np
from lxml import etree
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Directory of output figures
local_directory = "/home/alexandra/Documents/Data/Generic_expression_patterns/"
output_directory = "output/"
os.makedirs(output_directory, exist_ok=True)


# ## Function to plot

def make_figure_panel(filename, scale_x_input, scale_y_input, x_loc, y_loc):
    panel = sg.fromfile(filename)

    panel_size = (
        np.round(float(panel.root.attrib["width"][:-2]) * 1.33, 0),
        np.round(float(panel.root.attrib["height"][:-2]) * 1.33, 0),
    )

    scale_x = scale_x_input
    scale_y = scale_y_input

    print(f"original: {panel_size}")
    print(f"scaled:{(panel_size[0]*scale_x,panel_size[1]*scale_y)}")

    panel = panel.getroot()
    panel.scale_xy(x=scale_x, y=scale_y)
    panel.moveto(x_loc, y_loc)

    return panel


# ## Figure 1

# Create panels for figure 1
# Note panel 1A was made in google slides and exported from inkscape
panel_1a = make_figure_panel(
    "Fig-1A-SOPHIE workflow.svg", scale_x_input=4, scale_y_input=4, x_loc=30, y_loc=10
)
panel_1b = make_figure_panel(
    "../pseudomonas_analysis/pa_vae_clustering.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=500,
)
panel_1c = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=500,
    y_loc=500,
)

panel_1a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_1b_label = sg.TextElement(10, 500, "B", size=18, weight="bold", font="Verdana")
panel_1c_label = sg.TextElement(500, 500, "C", size=18, weight="bold", font="Verdana")

figure_1 = sg.SVGFigure("900", "900")
figure_1.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_1a,
        panel_1b,
        panel_1c,
        panel_1a_label,
        panel_1b_label,
        panel_1c_label,
    ]
)
display(SVG(figure_1.to_str()))

# save generated SVG files
figure_1.save("output/figure_1.svg")

# ## Figure 2

# Create panels for figure 2
panel_2a = make_figure_panel(
    "../human_cancer_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=10,
    y_loc=20,
)
panel_2b = make_figure_panel(
    "../human_general_analysis/gene_ranking_log2FoldChange.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=400,
    y_loc=20,
)
panel_2c = make_figure_panel(
    "../pseudomonas_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=800,
    y_loc=20,
)
panel_2d = make_figure_panel(
    "../compare_experiments/concordance_between_same_recount2_templates.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=30,
    y_loc=400,
)
panel_2e = make_figure_panel(
    "../compare_experiments/concordance_between_diff_recount2_templates.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=400,
    y_loc=400,
)

panel_2a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_2b_label = sg.TextElement(400, 20, "B", size=18, weight="bold", font="Verdana")
panel_2c_label = sg.TextElement(800, 20, "C", size=18, weight="bold", font="Verdana")
panel_2d_label = sg.TextElement(30, 400, "D", size=18, weight="bold", font="Verdana")
panel_2e_label = sg.TextElement(400, 400, "E", size=18, weight="bold", font="Verdana")

figure_2 = sg.SVGFigure("1200", "750")
figure_2.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_2a,
        panel_2b,
        panel_2c,
        panel_2d,
        panel_2e,
        panel_2a_label,
        panel_2b_label,
        panel_2c_label,
        panel_2d_label,
        panel_2e_label,
    ]
)
display(SVG(figure_2.to_str()))

# save generated SVG files
figure_2.save("output/figure_2.svg")

# ## Figure 3

# Create panels for figure 3
# Note: panel 3B was created in google slides and then exported from inkscape
panel_3a = make_figure_panel(
    "../human_cancer_analysis/pathway_ranking_padj.svg",
    scale_x_input=1.3,
    scale_y_input=1.3,
    x_loc=30,
    y_loc=10,
)
panel_3b = make_figure_panel(
    "Fig-4C-extendable workflow.svg",
    scale_x_input=4.5,
    scale_y_input=4.5,
    x_loc=600,
    y_loc=30,
)
panel_3c = make_figure_panel(
    "../other_enrichment_methods/enrichment_corr_plot.svg",
    scale_x_input=1.5,
    scale_y_input=1.5,
    x_loc=30,
    y_loc=400,
)


panel_3a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_3b_label = sg.TextElement(600, 30, "B", size=18, weight="bold", font="Verdana")
panel_3c_label = sg.TextElement(30, 400, "C", size=18, weight="bold", font="Verdana")

figure_3 = sg.SVGFigure("1600", "1000")
figure_3.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_3a,
        panel_3b,
        panel_3c,
        # panel_3d,
        panel_3a_label,
        panel_3b_label,
        panel_3c_label,
        # panel_3d_label,
    ]
)
display(SVG(figure_3.to_str()))

# save generated SVG files
figure_3.save("output/figure_3.svg")

# ## Figure 4

# Create panels for figure 4
panel_4a = make_figure_panel(
    "../LV_analysis/nonzero_LV_coverage.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=10,
)
panel_4b = make_figure_panel(
    "../LV_analysis/highweight_LV_coverage.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=350,
    y_loc=10,
)
panel_4c = make_figure_panel(
    "../LV_analysis/weight_dist_LV61.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=700,
    y_loc=10,
)
panel_4d = make_figure_panel(
    "../network_analysis/communities_fig.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=300,
)
panel_4e = make_figure_panel(
    "../network_analysis/centrality_figure.svg",
    scale_x_input=0.85,
    scale_y_input=0.85,
    x_loc=600,
    y_loc=300,
)

panel_4a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_4b_label = sg.TextElement(350, 20, "B", size=18, weight="bold", font="Verdana")
panel_4c_label = sg.TextElement(700, 20, "C", size=18, weight="bold", font="Verdana")
panel_4d_label = sg.TextElement(10, 300, "D", size=18, weight="bold", font="Verdana")
panel_4e_label = sg.TextElement(600, 300, "E", size=18, weight="bold", font="Verdana")

figure_4 = sg.SVGFigure("1050", "650")
figure_4.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_4a,
        panel_4b,
        panel_4c,
        panel_4d,
        panel_4e,
        panel_4a_label,
        panel_4b_label,
        panel_4c_label,
        panel_4d_label,
        panel_4e_label,
    ]
)
display(SVG(figure_4.to_str()))

# save generated SVG files
figure_4.save("output/figure_4.svg")

# ## Figure 5

# Create panels for figure 5
panel_5a = make_figure_panel(
    "Fig-5A_cbrAB_simpler_arg_model.svg",
    scale_x_input=2,
    scale_y_input=2,
    x_loc=30,
    y_loc=20,
)
panel_5b1 = make_figure_panel(
    "../pseudomonas_analysis/cbrB_volcano_zscore_highlight.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=450,
    y_loc=10,
)
panel_5b2 = make_figure_panel(
    "../pseudomonas_analysis/crc_volcano_zscore_highlight.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=850,
    y_loc=10,
)
panel_5c = make_figure_panel(
    "../pseudomonas_analysis/cbrB_crc_zscore_compare.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=450,
    y_loc=300,
)
panel_5d = make_figure_panel(
    os.path.join(local_directory, "WT-cbrB-EV-comp_10nM_Arg_crc.svg"),
    scale_x_input=2,
    scale_y_input=2,
    x_loc=1000,
    y_loc=300,
)

panel_5a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_5b_label = sg.TextElement(450, 20, "B", size=18, weight="bold", font="Verdana")
panel_5c_label = sg.TextElement(450, 300, "C", size=18, weight="bold", font="Verdana")
panel_5d_label = sg.TextElement(1000, 300, "D", size=18, weight="bold", font="Verdana")

figure_5 = sg.SVGFigure("1500", "700")
figure_5.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_5a,
        panel_5b1,
        panel_5b2,
        panel_5c,
        panel_5d,
        panel_5a_label,
        panel_5b_label,
        panel_5c_label,
        panel_5d_label,
    ]
)
display(SVG(figure_5.to_str()))

# save generated SVG files
figure_5.save("output/figure_5.svg")

# ## Supplement 1

# Create panels for Supplement 1
panel_S1a = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/array_expression_dist_gene_groups_highlight.svg",
    scale_x_input=0.9,
    scale_y_input=0.9,
    x_loc=30,
    y_loc=20,
)
panel_S1b = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/recount2_expression_dist_gene_groups_highlight.svg",
    scale_x_input=0.9,
    scale_y_input=0.9,
    x_loc=600,
    y_loc=10,
)
panel_S1c = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/violin_plot_mean_expression_RNAseq_only.svg",
    scale_x_input=0.6,
    scale_y_input=0.6,
    x_loc=30,
    y_loc=300,
)
panel_S1d = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/violin_plot_mean_expression_RNAseq_array.svg",
    scale_x_input=0.6,
    scale_y_input=0.6,
    x_loc=30,
    y_loc=1200,
)

panel_S1a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_S1b_label = sg.TextElement(600, 20, "B", size=18, weight="bold", font="Verdana")
panel_S1c_label = sg.TextElement(10, 300, "C", size=18, weight="bold", font="Verdana")
panel_S1d_label = sg.TextElement(30, 1200, "D", size=18, weight="bold", font="Verdana")

figure_S1 = sg.SVGFigure("1200", "2100")
figure_S1.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_S1a,
        panel_S1b,
        panel_S1c,
        panel_S1d,
        panel_S1a_label,
        panel_S1b_label,
        panel_S1c_label,
        panel_S1d_label,
    ]
)
display(SVG(figure_S1.to_str()))

# save generated SVG files
figure_S1.save("output/figure_S1.svg")

# ## Supplement 2

# Create panels for figure 2
panel_S2a = make_figure_panel(
    "../compare_experiments/cross_context_venn.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=10,
    y_loc=20,
)
panel_S2b = make_figure_panel(
    "../compare_experiments/cross_platform_venn.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=500,
    y_loc=20,
)

panel_S2a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_S2b_label = sg.TextElement(400, 20, "B", size=18, weight="bold", font="Verdana")

figure_S2 = sg.SVGFigure("950", "300")
figure_S2.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_S2a,
        panel_S2b,
        panel_S2a_label,
        panel_S2b_label,
    ]
)
display(SVG(figure_S2.to_str()))

# save generated SVG files
figure_S2.save("output/figure_S2.svg")

# ## Supplement 3

panel_Sa = make_figure_panel(
    "../other_enrichment_methods/enrichment_paired_plot_rnaseq.svg",
    scale_x_input=2,
    scale_y_input=2,
    x_loc=30,
    y_loc=20,
)

panel_Sa_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")

figure_S3 = sg.SVGFigure("1700", "1600")
figure_S3.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_Sa,
        panel_Sa_label,
    ]
)
display(SVG(figure_S3.to_str()))

# save generated SVG files
figure_S3.save("output/figure_S3.svg")

# ## Supplement 4

# Create panels for figure 2
panel_S4a = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC_5.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=10,
    y_loc=20,
)
panel_S4b = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC_10.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=400,
    y_loc=20,
)
panel_S4c = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC_25.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=800,
    y_loc=20,
)
panel_S4d = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC_50.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=30,
    y_loc=400,
)
panel_S4e = make_figure_panel(
    "../human_general_array_analysis/gene_ranking_logFC_100.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=400,
    y_loc=400,
)

panel_S4a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_S4b_label = sg.TextElement(400, 20, "B", size=18, weight="bold", font="Verdana")
panel_S4c_label = sg.TextElement(800, 20, "C", size=18, weight="bold", font="Verdana")
panel_S4d_label = sg.TextElement(30, 400, "D", size=18, weight="bold", font="Verdana")
panel_S4e_label = sg.TextElement(400, 400, "E", size=18, weight="bold", font="Verdana")

figure_S4 = sg.SVGFigure("1200", "750")
figure_S4.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_S4a,
        panel_S4b,
        panel_S4c,
        panel_S4d,
        panel_S4e,
        panel_S4a_label,
        panel_S4b_label,
        panel_S4c_label,
        panel_S4d_label,
        panel_S4e_label,
    ]
)
display(SVG(figure_S4.to_str()))

# save generated SVG files
figure_S4.save("output/figure_S4.svg")

# ## Output png version

# !inkscape --export-png=output/figure_1.png output/figure_1.svg
# !inkscape --export-png=output/figure_2.png output/figure_2.svg
# !inkscape --export-png=output/figure_3.png output/figure_3.svg
# !inkscape --export-png=output/figure_4.png output/figure_4.svg
# !inkscape --export-png=output/figure_5.png output/figure_5.svg
# !inkscape --export-png=output/figure_S1.png output/figure_S1.svg
# !inkscape --export-png=output/figure_S2.png output/figure_S2.svg
# !inkscape --export-png=output/figure_S3.png output/figure_S3.svg
# !inkscape --export-png=output/figure_S4.png output/figure_S4.svg
