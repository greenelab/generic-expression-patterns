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
panel_1 = make_figure_panel(
    "Fig-1 Overall workflow.svg", scale_x_input=4, scale_y_input=4, x_loc=30, y_loc=10
)

# +
# panel_1_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
# -

figure_1 = sg.SVGFigure("900", "500")
figure_1.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_1,
    ]
)
display(SVG(figure_1.to_str()))

# save generated SVG files
figure_1.save("output/figure_1.svg")

# ## Figure 2

# Create panels for figure 2
panel_2a = make_figure_panel(
    "Fig-2A-genes workflow.svg", scale_x_input=4, scale_y_input=4, x_loc=30, y_loc=10
)
panel_2b = make_figure_panel(
    "Fig-2B-model validation.svg",
    scale_x_input=2,
    scale_y_input=2,
    x_loc=30,
    y_loc=550,
)
panel_2c = make_figure_panel(
    "Fig-2C-model power.svg",
    scale_x_input=2,
    scale_y_input=2,
    x_loc=400,
    y_loc=550,
)

panel_2a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_2b_label = sg.TextElement(10, 550, "B", size=18, weight="bold", font="Verdana")
panel_2c_label = sg.TextElement(400, 550, "C", size=18, weight="bold", font="Verdana")

figure_2 = sg.SVGFigure("1000", "1000")
figure_2.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_2a,
        panel_2b,
        panel_2c,
        panel_2a_label,
        panel_2b_label,
        panel_2c_label,
    ]
)
display(SVG(figure_2.to_str()))

# save generated SVG files
figure_2.save("output/figure_2.svg")

# ## Figure 3

# Create panels for figure 2
panel_3a = make_figure_panel(
    "../human_cancer_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=10,
    y_loc=20,
)
panel_3b = make_figure_panel(
    "../human_general_analysis/gene_ranking_log2FoldChange.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=400,
    y_loc=20,
)
panel_3c = make_figure_panel(
    "../pseudomonas_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=800,
    y_loc=20,
)
panel_3d = make_figure_panel(
    "../compare_experiments/concordance_between_same_recount2_templates.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=30,
    y_loc=400,
)
panel_3e = make_figure_panel(
    "../compare_experiments/concordance_between_diff_recount2_templates.svg",
    scale_x_input=0.77,
    scale_y_input=0.77,
    x_loc=400,
    y_loc=400,
)

panel_3a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_3b_label = sg.TextElement(400, 20, "B", size=18, weight="bold", font="Verdana")
panel_3c_label = sg.TextElement(800, 20, "C", size=18, weight="bold", font="Verdana")
panel_3d_label = sg.TextElement(30, 400, "D", size=18, weight="bold", font="Verdana")
panel_3e_label = sg.TextElement(400, 400, "E", size=18, weight="bold", font="Verdana")

figure_3 = sg.SVGFigure("1200", "800")
figure_3.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_3a,
        panel_3b,
        panel_3c,
        panel_3d,
        panel_3e,
        panel_3a_label,
        panel_3b_label,
        panel_3c_label,
        panel_3d_label,
        panel_3e_label,
    ]
)
display(SVG(figure_3.to_str()))

# save generated SVG files
figure_3.save("output/figure_3.svg")

# ## Figure 4

# Create panels for figure 3
panel_4a = make_figure_panel(
    "Fig-4A-pathway workflow .svg", scale_x_input=4, scale_y_input=4, x_loc=30, y_loc=10
)
panel_4b = make_figure_panel(
    "../human_cancer_analysis/pathway_ranking_padj.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=900,
    y_loc=30,
)
panel_4c = make_figure_panel(
    "Fig-4C-extendable workflow.svg",
    scale_x_input=3,
    scale_y_input=3,
    x_loc=30,
    y_loc=600,
)
panel_4d = make_figure_panel(
    "../other_enrichment_methods/enrichment_corr_plot.svg",
    scale_x_input=1.5,
    scale_y_input=1.5,
    x_loc=750,
    y_loc=600,
)


panel_4a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_4b_label = sg.TextElement(900, 20, "B", size=18, weight="bold", font="Verdana")
panel_4c_label = sg.TextElement(10, 600, "C", size=18, weight="bold", font="Verdana")
panel_4d_label = sg.TextElement(700, 600, "D", size=18, weight="bold", font="Verdana")

figure_4 = sg.SVGFigure("1500", "1200")
figure_4.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_4a,
        panel_4b,
        panel_4c,
        panel_4d,
        panel_4a_label,
        panel_4b_label,
        panel_4c_label,
        panel_4d_label,
    ]
)
display(SVG(figure_4.to_str()))

# save generated SVG files
figure_4.save("output/figure_4.svg")

# ## Figure 5

# Create panels for figure 4
panel_5a = make_figure_panel(
    "../LV_analysis/nonzero_LV_coverage.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=10,
)
panel_5b = make_figure_panel(
    "../LV_analysis/highweight_LV_coverage.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=350,
    y_loc=10,
)
panel_5c = make_figure_panel(
    "../LV_analysis/weight_dist_LV61.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=700,
    y_loc=10,
)
panel_5d = make_figure_panel(
    "../network_analysis/communities_fig.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=300,
)
panel_5e = make_figure_panel(
    "../network_analysis/centrality_figure.svg",
    scale_x_input=0.85,
    scale_y_input=0.85,
    x_loc=600,
    y_loc=300,
)

panel_5a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_5b_label = sg.TextElement(350, 20, "B", size=18, weight="bold", font="Verdana")
panel_5c_label = sg.TextElement(700, 20, "C", size=18, weight="bold", font="Verdana")
panel_5d_label = sg.TextElement(10, 300, "D", size=18, weight="bold", font="Verdana")
panel_5e_label = sg.TextElement(600, 300, "E", size=18, weight="bold", font="Verdana")

figure_5 = sg.SVGFigure("1500", "800")
figure_5.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_5a,
        panel_5b,
        panel_5c,
        panel_5d,
        panel_5e,
        panel_5a_label,
        panel_5b_label,
        panel_5c_label,
        panel_5d_label,
        panel_5e_label,
    ]
)
display(SVG(figure_5.to_str()))

# save generated SVG files
figure_5.save("output/figure_5.svg")

# ## Figure 6

# Create panels for figure 5
panel_6a = make_figure_panel(
    os.path.join(local_directory, "cbrAB_simpler_arg_model.svg"),
    scale_x_input=2,
    scale_y_input=2,
    x_loc=30,
    y_loc=20,
)
panel_6b = make_figure_panel(
    os.path.join(local_directory, "template_zscore_volcano_ArgR_E-GEOD-33245.svg"),
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=350,
    y_loc=10,
)
panel_6c = make_figure_panel(
    os.path.join(local_directory, "template_traditional_volcano_ArgR_E-GEOD-33245.svg"),
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=750,
    y_loc=10,
)
panel_6d = make_figure_panel(
    "../pseudomonas_analysis/cbrB_crc_zscore_compare.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=350,
    y_loc=300,
)
panel_6e = make_figure_panel(
    os.path.join(local_directory, "WT-cbrB-EV-comp_10nM_Arg_crc.svg"),
    scale_x_input=2,
    scale_y_input=2,
    x_loc=750,
    y_loc=300,
)

panel_6a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_6b_label = sg.TextElement(350, 20, "B", size=18, weight="bold", font="Verdana")
panel_6c_label = sg.TextElement(750, 20, "C", size=18, weight="bold", font="Verdana")
panel_6d_label = sg.TextElement(350, 300, "D", size=18, weight="bold", font="Verdana")
panel_6e_label = sg.TextElement(750, 300, "E", size=18, weight="bold", font="Verdana")

figure_6 = sg.SVGFigure("1200", "800")
figure_6.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_6a,
        panel_6b,
        panel_6c,
        panel_6d,
        panel_6e,
        panel_6a_label,
        panel_6b_label,
        panel_6c_label,
        panel_6d_label,
        panel_6e_label,
    ]
)
display(SVG(figure_6.to_str()))

# save generated SVG files
figure_6.save("output/figure_6.svg")

# ## Supplement 1

# Create panels for Supplement 1
panel_S1a = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/array_expression_dist_gene_groups_highlight.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=30,
    y_loc=20,
)
panel_S1b = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/recount2_expression_dist_gene_groups_highlight.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=600,
    y_loc=10,
)
panel_S1c = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/violin_plot_mean_expression_RNAseq_only.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=300,
)
panel_S1d = make_figure_panel(
    "../explore_RNAseq_only_generic_genes/violin_plot_mean_expression_RNAseq_array.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=1800,
)

panel_S1a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_S1b_label = sg.TextElement(600, 20, "B", size=18, weight="bold", font="Verdana")
panel_S1c_label = sg.TextElement(10, 300, "C", size=18, weight="bold", font="Verdana")
panel_S1d_label = sg.TextElement(30, 1800, "D", size=18, weight="bold", font="Verdana")

figure_S1 = sg.SVGFigure("1600", "3100")
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

panel_Sa = make_figure_panel(
    "../other_enrichment_methods/enrichment_paired_plot_rnaseq.svg",
    scale_x_input=2,
    scale_y_input=2,
    x_loc=30,
    y_loc=20,
)

panel_Sa_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")

figure_S2 = sg.SVGFigure("1800", "1800")
figure_S2.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_Sa,
        panel_Sa_label,
    ]
)
display(SVG(figure_S2.to_str()))

# save generated SVG files
figure_S2.save("output/figure_S2.svg")

# ## Supplement 3

# Create panels for Supplement 2
panel_S3a = make_figure_panel(
    "../human_general_analysis/gene_ranking_log2FoldChange.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=20,
)
panel_S3b = make_figure_panel(
    "../human_cancer_analysis/gene_ranking_logFC.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=400,
    y_loc=10,
)
panel_S3c = make_figure_panel(
    "../human_general_analysis/pathway_ranking_padj.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=30,
    y_loc=400,
)
panel_S3d = make_figure_panel(
    "../human_cancer_analysis/pathway_ranking_padj.svg",
    scale_x_input=0.8,
    scale_y_input=0.8,
    x_loc=400,
    y_loc=400,
)

panel_S3a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_S3b_label = sg.TextElement(400, 20, "B", size=18, weight="bold", font="Verdana")
panel_S3c_label = sg.TextElement(10, 400, "C", size=18, weight="bold", font="Verdana")
panel_S3d_label = sg.TextElement(400, 400, "D", size=18, weight="bold", font="Verdana")

figure_S3 = sg.SVGFigure("800", "800")
figure_S3.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_S3a,
        panel_S3b,
        panel_S3c,
        panel_S3d,
        panel_S3a_label,
        panel_S3b_label,
        panel_S3c_label,
        panel_S3d_label,
    ]
)
display(SVG(figure_S3.to_str()))

# save generated SVG files
figure_S3.save("output/figure_S3.svg")

# ## Output png version

# !inkscape --export-png=output/figure_1.png output/figure_1.svg
# !inkscape --export-png=output/figure_2.png output/figure_2.svg
# !inkscape --export-png=output/figure_3.png output/figure_3.svg
# !inkscape --export-png=output/figure_4.png output/figure_4.svg
# !inkscape --export-png=output/figure_5.png output/figure_5.svg
# !inkscape --export-png=output/figure_6.png output/figure_6.svg
# !inkscape --export-png=output/figure_S1.png output/figure_S1.svg
# !inkscape --export-png=output/figure_S2.png output/figure_S2.svg
# !inkscape --export-png=output/figure_S3.png output/figure_S3.svg
