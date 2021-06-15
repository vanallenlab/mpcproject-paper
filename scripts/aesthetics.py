"""Shared settings for controlling figure aesthetics. Taken with permission from Meng Xiao He"""
import functools

import palettable
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib import rcParams


custom_rcParams = {
    'figure.figsize': (8, 3),
    'font.family': 'Gudea',
    'font.size': 12,
    'font.weight': 'regular',
    'axes.labelsize': 13,
    'axes.formatter.useoffset': False,
    'axes.formatter.limits': (-4, 4),
    'axes.titlesize': 14,
    'legend.edgecolor': 'none',
    'legend.fancybox': False,
    'legend.fontsize': 13,
    'legend.frameon': False,
    'legend.framealpha': 0,
    'legend.facecolor': 'none',
    'legend.loc': 'center left',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
}

rcParams.update(custom_rcParams)
sns.set_context('paper', rc=custom_rcParams)
# set legend to be outside of plot area by default
plt.legend = functools.partial(plt.legend, bbox_to_anchor=(1, 0.5))


ASPECT_RATIO = 4 / 3

paper_rcParams = {
    **custom_rcParams,
    'axes.prop_cycle': cycler(color=palettable.cartocolors.qualitative.Bold_10.mpl_colors),
    'figure.figsize': (6, 6),
    'figure.dpi': 96,
    'font.family': 'Arial',
    'legend.fontsize': 12,
    'savefig.dpi': 300
}

def activate_paper_rcParams():
    rcParams.update(paper_rcParams)
    sns.set_context('paper', rc=paper_rcParams)
