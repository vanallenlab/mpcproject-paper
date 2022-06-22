import numpy as np
import seaborn as sns
from matplotlib_venn import venn3, venn3_circles
import math

# monkey patch stripplot to be deterministic
original_stripplot = sns.stripplot
def stripplot(*args, **kwargs):
    random_state = np.random.get_state()
    np.random.seed(42)
    p = original_stripplot(*args, **kwargs)
    np.random.set_state(random_state)
    return p
sns.stripplot = stripplot

# a minor formatting function
def format_func(value, tick_number):
    if value == 0:
        return "0"
    elif value == 1:
        return r"$\mathregular{10^{1}}$"
    elif value == 2:
        return r"$\mathregular{10^{2}}$"
    elif value == 3:
        return r"$\mathregular{10^{3}}$"
    elif value == 4:
        return r"$\mathregular{10^{4}}$"

def venn_diagram(a, b, c, ax, labels=['A', 'B', 'C'], circles = True):

    a = set(a)
    b = set(b)
    c = set(c)

    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)

    only_a_b = len(a & b - c)
    only_a_c = len(a & c - b)
    only_b_c = len(b & c - a)

    a_b_c = len(a & b & c)
    
    if circles:
        v = venn3_circles(subsets=(only_a, only_b, only_a_b, only_c, only_a_c, only_b_c, a_b_c), 
                          ax = ax)
    else:
        v = venn3(subsets=(only_a, only_b, only_a_b, only_c, only_a_c, only_b_c, a_b_c), 
                  ax=ax, set_labels=labels)
    return v
    
def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

# A custom function to interpolate dates
def date2plotcoords(date1_tuple, date2_tuple):
    '''Uses 2 "anchors" to define the slope and intercept of a date to plot
    coordinate conversion.'''
    
    a = (date1_tuple[1] - date2_tuple[1])/(date1_tuple[0] - date2_tuple[0])
    b = date1_tuple[1] - a*date1_tuple[0]
    
    return lambda date: (date - b)/a

def clean_gistic_file(gistic_file, amp_or_del):
    '''Takes in a GISTIC amp or del text file and cleans it for CBIO release'''
    
    # import file and transpose, as GISTIC has row names
    gistic_file = gistic_file.transpose().iloc[1:, :].copy()

    # drop rows with all NAs (extra rows)
    gistic_file.dropna(how = 'all', inplace = True)
    
    # the first four columns contain info about regions
    region_info = gistic_file[gistic_file.columns[:4]].copy()
    region_info.columns = ['cytoband', 'q_value', 'residual_q_value', 'cna_region']

    # while the 5th column onwards contains info about genes
    genes = gistic_file.loc[:, gistic_file.columns[4:]].copy()

    # collapse the genes into one comma separated column, then add back to region_info
    gene_names = genes.apply(lambda x: ','.join(x.dropna()), axis = 1)
    gene_lengths = genes.apply(lambda x: len(x.dropna()), axis = 1)
    region_info.loc[:, 'genes_in_region'] = gene_names.copy()
    region_info.loc[:, 'n_genes_in_region'] = gene_lengths.copy()

    # add necessary columns - chromosome, peak start/stop, n_genes, and amp
    region_info.loc[:, 'chromosome'] = region_info['cna_region'].apply(lambda x: x.split(':')[0][3:])
    coordinates = region_info['cna_region'].str.split(':').str[1]
    region_info[['peak_start', 'peak_end']] = coordinates.str.split('-', expand = True).copy(deep = True)
    region_info.loc[:, 'amp'] = 1 if amp_or_del == 'amp' else 0

    # reorder columns and export
    out_df = region_info[['chromosome', 'peak_start', 'peak_end', 'n_genes_in_region', 'genes_in_region', 'cytoband', 'q_value', 'residual_q_value', 'amp']]

    return out_df