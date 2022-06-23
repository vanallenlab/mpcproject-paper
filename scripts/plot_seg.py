import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable

def transform_to_int(object):
    if isinstance(object, (int, float)): 
        return object
    elif isinstance(object, str):
        if object.isdigit():
            return int(object)
        else:
            return object
    else:
        raise ValueError('Unrecognized chromosome {} of type {}'.format(object, type(object)))


class SegObject:

    '''A class designed to handle segmentation files for visualization of copy number'''

    def __init__(self, sampleid, seg_df=None, seg_path=None, seg_type=None, start_name=None, end_name=None,
                 chrom_name=None, values=None, chrom_df=None, **kwargs):

        '''Params:
           -------
           seg_path: path to seg file (either this or seg_df must be defined)
           seg_type: 'facets', 'absolute', 'gatk'
           values: list of column names to plot (length 1 for total, 2 for allelic), must be provided if seg_type is not passed
           start_name, end_name, chrom_name: can be used to override defaults, must be provided if seg_type is not passed
           kwargs: kwargs to pass to pd.read_csv'''

        self.sampleid = sampleid

        # attributes of the seg file
        if seg_type:
            if seg_type == 'facets':
                start_name, end_name, chrom_name = 'start', 'end', 'chrom'
                values = ['tcn.em', 'lcn.em']
            elif seg_type == 'absolute':
                start_name, end_name, chrom_name = 'Start.bp', 'End.bp', 'Chromosome'
                values = ['modal.a1', 'modal.a2']
            elif seg_type == 'gatk':
                start_name, end_name, chrom_name = 'Start', 'End', 'Chromosome'
                values = ['Segment_Mean']
            else:
                raise ValueError('seg type not recognized. Valid options are facets, absolute, or gatk')

        else:
            if any([value is None for value in [start_name, end_name, chrom_name, values]]):
                raise ValueError('start_name, end_name, chrom_name, and values must be defined if no seg_type passed')

        self.start_name = start_name
        self.end_name = end_name
        self.chrom_name = chrom_name
        self.values = values

        # the seg data
        if seg_df is None:
            if seg_path is None:
                raise ValueError('seg_path must be provided if seg_df is not')
            seg_df = pd.read_csv(seg_path, **kwargs)

        if not set([start_name, end_name, chrom_name]).issubset(set(seg_df.columns)):
            diff = set([start_name, end_name, chrom_name]) - set(seg_df.columns)
            raise ValueError('{} not in seg columns'.format(diff))

        self.seg_data = seg_df

        # set the reference chromosome locations
        if chrom_df is None:
            # defaults taken from the broad's reference genome hg19
            chroms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']
            lengths = [249250621.0, 243199373.0, 198022430.0, 191154276.0, 180915260.0, 171115067.0, 159138663.0, 146364022.0, 141213431.0, 135534747.0, 
                    135006516.0, 133851895.0, 115169878.0, 107349540.0, 102531392.0, 90354753.0, 81195210.0, 78077248.0, 59128983.0, 63025520.0, 48129895.0, 
                    51304566.0, 155270560.0, 59373566.0]
            centromeres = [125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000, 40200000, 53700000, 35800000, 17900000, 17600000, 
                           19000000, 36600000, 24000000, 17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 12500000]
            chrom_df = pd.DataFrame(list(zip(chroms, lengths)), columns =['chrom', 'length'])
            chrom_df['centromere'] = centromeres

        # define the cumulative sum
        chrom_df['cum_sum'] = [0] + list(np.cumsum(chrom_df['length'][:-1]))
        chrom_dict = dict(zip(chrom_df['chrom'], chrom_df['cum_sum']))

        self.chrom_df = chrom_df
        self.chrom_dict = chrom_dict

    def plot_seg_data(self, ax, ylim=None, **kwargs):
        '''Plots the data associated with the segmentation file.'''

        # plot the "background"
        colors = ['white', 'grey']
        start = 0
        xticks = []
        for index, row in self.chrom_df.iterrows():
            end = start + row['length']
            centromere = start + row['centromere']
            middle = (start + end)/2
            ax.axvspan(start, end, color = colors[index % 2], alpha = 0.25)
            ax.axvline(x = centromere, linestyle = 'dotted', color = 'black', alpha = 0.75)
            start += row['length']
            xticks.append((row['chrom'], middle))

        # define colors for the plot
        if len(self.values) == 2:
            colors = ['red', 'blue']
        elif len(self.values) == 1:
            colors = ['black']
        else:
            raise ValueError('Values must be length 1 (total cn) or 2 (allelic cn)')

        # now plot the data
        for index, seg in self.seg_data.iterrows():
            start, end, chrom = seg[self.start_name], seg[self.end_name], seg[self.chrom_name]
            copy_nums = [seg[val] for val in self.values]
            if not all([isinstance(cn, (float, int)) for cn in copy_nums]):
                raise ValueError('At least one of {} is not a numeric copy number'.format(copy_nums))

            if chrom == 23: chrom = 'X'
            if chrom == 24: chrom = 'Y'

            # plot total or allelic copy number
            chrom_length = self.chrom_dict[chrom]
            for i, cn in enumerate(copy_nums):
                ax.plot([chrom_length + start, chrom_length + end], [cn, cn], color = colors[i], alpha = 0.5, **kwargs)

        ax.set_xlim(0, end)
        if not ylim is None:
            ax.set_ylim(ylim[0], ylim[1])
        ax.set_xticks([loc for _, loc in xticks])
        ax.set_xticklabels([label for label, _ in xticks])
        ax.tick_params(axis = 'x', length = 0)

        ax.set_title(self.sampleid)
        return ax

    def visualize_cnas(self, ax, total_cn_col, color_dict=None, orientation='vertical'):
        '''Visualizes the total copy nummber in a nice summarized plot.
        Does not correctly handle WGD samples, as it visualizes raw CN. Does not handle CN > 4.

        ax: axis object to plot on
        total_cn_col: column name of the total copy number, used to shade
        color_dict: mapping of copy number to color. Defaults to palettable (see below)
        orientation: orientation of the plot.
        '''

        if color_dict is None:

            # not super satisfied with the 4 CN
            balance_6 = palettable.cmocean.diverging.Balance_6.mpl_colors
            color_dict = {0: {'facecolor': (0,0,0), 'alpha': 0.5}, 1: {'facecolor': balance_6[1], 'alpha': 0.5},
                         2: {'facecolor': (1,1,1)}, 3: {'facecolor': balance_6[-2], 'alpha': 0.5}}

            # an ugly hack to introduce higher amplifications
            for cn in range(4, 100):
                color_dict[cn] = {'facecolor': balance_6[-2]}

        self.color_dict = color_dict

        for index, seg in self.seg_data.iterrows():

            start, end, chrom = seg[self.start_name], seg[self.end_name], seg[self.chrom_name]
            total_cn = seg[total_cn_col]

            if chrom == 23: chrom = 'X'
            if chrom == 24: chrom = 'Y'

            # plot total copy number
            chrom_start = self.chrom_dict[chrom]
            if orientation == 'vertical':
                ax.axhspan(ymin = chrom_start + start, ymax = chrom_start + end, 
                        xmin = 0, xmax = 1, **color_dict[total_cn])

            else:
                ax.axvspan(xmin = chrom_start + start, xmax = chrom_start + end, 
                        ymin = 0, ymax = 1, **color_dict[total_cn])

        # plot chromosome breaks
        labels, locs = [], []
        for index, row in self.chrom_df.iterrows():
            chrom = row['chrom']
            if not chrom in ['X', 'Y']:
                labels.append(chrom)
                locs.append(row['cum_sum'] + row['length']/2)

                # no need for division at chrom 1
                if chrom == 1: continue
                if orientation == 'vertical':
                    ax.plot([0, 1], [row['cum_sum'], row['cum_sum']], color = 'black')
                else:
                    ax.plot([row['cum_sum'], row['cum_sum']], [0, 1], color = 'black', linestyle = 'dotted')

        if orientation == 'vertical':
            ax.set_yticks(locs)
            ax.set_yticklabels(labels)
            ax.set_xlim([0, 1])
            ax.set_ylim([0, self.chrom_dict['X']])
            ax.get_xaxis().set_visible(False)
            ax.tick_params(axis = 'y', length=0)
            ax.set_title(self.sampleid)

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False) 
            ax.spines['bottom'].set_visible(False)
        else:
            ax.set_xticks(locs)
            ax.set_xticklabels(labels)
            ax.set_ylim([0, 1])
            ax.set_xlim([0, self.chrom_dict['X']])
            ax.set_ylabel(self.sampleid)
            ax.set_yticklabels([])
            ax.tick_params(axis = 'x', length=0)
            ax.tick_params(axis = 'y', length=0)

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False) 
            ax.spines['left'].set_visible(False)

        return ax
