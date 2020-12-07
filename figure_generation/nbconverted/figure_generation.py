#!/usr/bin/env python
# coding: utf-8

# # Figure generation

# In[1]:


from IPython.display import Image, display, SVG
import svgutils.transform as sg
import numpy as np
from lxml import etree
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# In[2]:


# Directory of output figures
output_directory = "output/"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)


# ## Function to plot

# In[3]:


def make_figure_panel(filename, scale_x_input, scale_y_input, x_loc, y_loc):
    panel = (
        sg.fromfile(filename)
    )
    
    panel_size = (
        np.round(float(panel.root.attrib['width'][:-2])*1.33, 0), 
        np.round(float(panel.root.attrib['height'][:-2])*1.33, 0)
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

# In[4]:


# Create panels for figure 1
panel_1a = make_figure_panel(
    "fig1A.svg",
    2.5, 
    2.5,
    30,
    10
)
panel_1b = make_figure_panel(
    "../human_general_analysis/logs/NN_2500_30/tybalt_2layer_30latent_hist.svg",
    0.8,
    0.8,
    600,
    10
)
panel_1c = make_figure_panel(
    "fig1C.svg",
    4.2,
    4.2,
    30,
    250
)
## TO DO
# Add panel showing volcano plots of example simulated experiments
# Annotate volcano with name of DEGs


# In[5]:


panel_1a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_1b_label = sg.TextElement(580, 20, "B", size=18, weight="bold", font="Verdana")
panel_1c_label = sg.TextElement(10, 260, "C", size=18, weight="bold", font="Verdana")


# In[6]:


figure_1 = sg.SVGFigure("1000", "1000")
figure_1.append([
    etree.Element("rect", {"width":"100%", "height":"100%", "fill":"white"}),
    panel_1a, 
    panel_1b,
    panel_1c,
    panel_1a_label,
    panel_1b_label,
    panel_1c_label
])
display(SVG(figure_1.to_str()))


# In[7]:


# save generated SVG files
figure_1.save("output/figure_1.svg")


# ## Figure 2

# In[8]:


# Create panels for figure 2
panel_2a = make_figure_panel(
    "fig2A.svg",
    2,
    2,
    30,
    10
)
panel_2b = make_figure_panel(
    "../human_cancer_analysis/gene_ranking_logFC.svg",
    0.6,
    0.6,
    470,
    30
)
panel_2c = make_figure_panel(
    "../human_general_analysis/gene_ranking_log2FoldChange.svg",
    0.6,
    0.6,
    700,
    30
)
panel_2d = make_figure_panel(
    "fig2D.svg",
    2,
    2,
    30,
    300
)
panel_2e = make_figure_panel(
    "../human_cancer_analysis/pathway_ranking_padj.svg",
    0.6,
    0.6,
    470,
    300
)


# In[9]:


panel_2a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_2b_label = sg.TextElement(450, 20, "B", size=18, weight="bold", font="Verdana")
panel_2c_label = sg.TextElement(680, 20, "C", size=18, weight="bold", font="Verdana")
panel_2d_label = sg.TextElement(10, 310, "D", size=18, weight="bold", font="Verdana")
panel_2e_label = sg.TextElement(450, 310, "E", size=18, weight="bold", font="Verdana")


# In[10]:


figure_2 = sg.SVGFigure("1000", "800")
figure_2.append([
    etree.Element("rect", {"width":"100%", "height":"100%", "fill":"white"}),
    panel_2a, 
    panel_2b,
    panel_2c,
    panel_2d,
    panel_2e,
    panel_2a_label,
    panel_2b_label,
    panel_2c_label,
    panel_2d_label,
    panel_2e_label
])
display(SVG(figure_2.to_str()))


# In[11]:


# save generated SVG files
figure_2.save("output/figure_2.svg")


# ## Figure 3

# In[12]:


# Create panels for figure 3
panel_3a = make_figure_panel(
    "../multiplier_analysis/nonzero_LV_coverage.svg",
    1,
    1,
    30,
    10
)
panel_3b = make_figure_panel(
    "../multiplier_analysis/highweight_LV_coverage.svg",
    1,
    1,
    500,
    10
)
## TO DO
# Add panel showing proportion of generic genes per LV


# In[13]:


panel_3a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_3b_label = sg.TextElement(480, 20, "B", size=18, weight="bold", font="Verdana")


# In[14]:


figure_3 = sg.SVGFigure("1000", "500")
figure_3.append([
    etree.Element("rect", {"width":"100%", "height":"100%", "fill":"white"}),
    panel_3a, 
    panel_3b,
    panel_3a_label,
    panel_3b_label,
])
display(SVG(figure_3.to_str()))


# In[15]:


# save generated SVG files
figure_3.save("output/figure_3.svg")

