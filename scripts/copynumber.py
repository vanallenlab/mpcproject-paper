import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable
import seaborn as sns
import warnings


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


def shift_loc(chrom_dict, chrom_list, loc_list):
    '''A function to shift positions to match a global profile'''

    offsets = np.array([chrom_dict[chrom] for chrom in chrom_list])
    return offsets + np.array(loc_list)


class CopyNumber:

    '''A class designed to handle data for visualization of copy number'''

    def __init__(self, sampleid, chrom_df=None):

        '''
        Params:
        -------
        sampleid: name of sample in question
        chrom_df: dataframe detailing the name, length, and centromere location of
            chromosomes. Defaults to the Broad's hg19.

        Attributes:
        -----------
        sampleid: the sample id of the sample being visualized
        chrom_df: dataframe containing chromosome info
        segments: pandas dataframe, containing segmentation data
        segment_values: list, the columns that will be plotted for segments
        coverage: pandas dataframe, containing coverage data
        '''

        self.sampleid = sampleid

        # handle default chrom_df
        if chrom_df is None:
            chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
            lengths = [249250621.0, 243199373.0, 198022430.0, 191154276.0, 180915260.0, 171115067.0, 159138663.0, 146364022.0, 141213431.0, 135534747.0, 
                    135006516.0, 133851895.0, 115169878.0, 107349540.0, 102531392.0, 90354753.0, 81195210.0, 78077248.0, 59128983.0, 63025520.0, 48129895.0, 
                    51304566.0, 155270560.0, 59373566.0]
            centromeres = [125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000, 40200000, 53700000, 35800000, 17900000, 17600000, 
                           19000000, 36600000, 24000000, 17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 12500000]
            chrom_df = pd.DataFrame(list(zip(chroms, lengths, centromeres)),
                                    columns =['chrom', 'length', 'centromere'])

            chrom_cum_sum = [0] + list(np.cumsum(chrom_df['length'][:-1]))
            chrom_dict = dict(zip(chrom_df['chrom'], chrom_cum_sum))
            self._chrom_dict = chrom_dict
        self.chrom_df = chrom_df

        self.segments = None
        self.segment_values = None
        self.coverage = None
        self.coverage_values = None

    def set_segments(self, seg_df, seg_values):

        '''Sets the segment and segment values attributes

        Params:
        -------
        seg_df: pandas dataframe
            DataFrame containing the segmentation data. Required columns are
            start, end, chrom, and the columns passed as segment_values

        seg_values: list
            List containing column names for the value that should be plotted.
            If more than one value is provided, multiple segments may be plotted.
        '''

        required_cols = ['start', 'end', 'chrom'] + list(seg_values)
        missing_cols = set(required_cols) - set(seg_df.columns)
        if missing_cols:
            raise ValueError('Missing columns: {}'.format(missing_cols))

        self.segments = seg_df
        self.segment_values = seg_values

        return self

    def set_coverage(self, cov_df, cov_values, downsample = 1):

        '''Sets the segment and segment values attributes

        Params:
        -------
        cov_df: pandas dataframe
            DataFrame containing the coverage data. Required columns are
            chrom, loc, and the columns passed as cov_values.

        cov_values: list
            List containing column names for the value that should be plotted.
            If more than one value is provided, multiple scatterplots may be plotted.
            
        downsample: int, default 1
            Whether or not to downsample coverage. Every n position will be included.
            e.g. if downsample is 2, every other position is included, etc.
        '''

        required_cols = ['chrom', 'loc'] + list(cov_values)
        missing_cols = set(required_cols) - set(cov_df.columns)
        if missing_cols:
            raise ValueError('Missing columns: {}'.format(missing_cols))
        
        cov_df = cov_df.iloc[::downsample, :]
        self.coverage = cov_df
        self.coverage_values = cov_values

        return self

                
    def visualize_profile(self, ax, coverage=True, segments=True, window=None, cov_colors=None,
                          ylim=None, seg_colors=None, plot_kwargs=None, scatter_kwargs=None):
        '''Visualizes the coverage and segment copy number on a plot.

        Params:
        -------
        ax: axis object
            Axis object on which to put the plot
        coverage: bool, default True
            Whether coverage should be plotted
        segments: bool, default True
            Whether segments should be plotted
        window: list, default None
            If provided, a window on which to "zoom" in.
            Format is (chr, start, end)
        cov_colors, seg_colors: list
            List of colors to be given for the coverage and segmentation
            values.
        ylim: list
            y limits for the plot (defaults to matplotlib's calculated values)
        plot_kwargs: dict, default None
            dict of kwargs to be passed to ax.plot in plotting
            segments
        scatter_kwargs: dict, default None
            dict of kwargs to be passed to sns.scatterplot in
            plotting coverage
        '''

        # make copies
        if coverage:
            cov_df = self.coverage.copy()
        if segments:
            seg_df = self.segments.copy()

        # handle defaults
        if cov_colors is None and coverage:
            default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            cov_colors = default_colors[:len(self.coverage_values)]

        if seg_colors is None and segments:
            seg_colors = ['black']*len(self.segment_values)

        if scatter_kwargs is None:
            scatter_kwargs = {'edgecolor': None, 's': 5}

        if plot_kwargs is None:
            plot_kwargs = {}

        # if there is no window, we plot the background
        if window is None:
            colors = ['white', 'grey']
            start = 0
            xticks = []
            for index, row in self.chrom_df.iterrows():
                end = start + row['length']
                centromere = start + row['centromere']
                middle = (start + end)/2
                ax.axvspan(start, end, color=colors[index % 2], alpha=0.25)
                ax.axvline(x=centromere, linestyle='dotted', color='black', alpha=0.75)
                start += row['length']
                xticks.append((row['chrom'], middle))

        # plot the coverage data, if it is exists
        if coverage:
            if window is None:
                # shift all coverage points to match global profile
                cov_df['loc'] = shift_loc(self._chrom_dict, cov_df['chrom'], cov_df['loc'])

            else:
                # find coverage in window
                window_chrom, window_start, window_end = window
                cov_df = cov_df[(cov_df['chrom'] == window_chrom) &
                                (cov_df['loc'] > window_start) &
                                (cov_df['loc'] < window_end)]
                if cov_df.empty:
                    warnings.warn("No coverage found in window {}:{}-{}".format(window_chrom, window_start, window_end))

            # plot all the values provided
            for i, value in enumerate(self.coverage_values):
                sns.scatterplot(data=cov_df, x='loc', y=value, ax=ax, color=cov_colors[i],
                                **scatter_kwargs)

        # plot the segment data, if it exists
        if segments:
            if window is None:

                # shift all segments to match global profile
                seg_df['start'] = shift_loc(self._chrom_dict, seg_df['chrom'], seg_df['start'])
                seg_df['end'] = shift_loc(self._chrom_dict, seg_df['chrom'], seg_df['end'])

            else:
                # find segments in window
                window_chrom, window_start, window_end = window
                seg_df = seg_df[seg_df['chrom'] == window_chrom]
                seg_df = seg_df[((seg_df['start'] <= window_start) & (seg_df['end'] >= window_start)) | 
                                ((seg_df['start'] >= window_start) & (seg_df['end'] <= window_end)) | 
                                ((seg_df['start'] <= window_end) & (seg_df['end'] >= window_end))]
                if seg_df.empty:
                    warnings.warn("No segments found in window {}:{}-{}".format(window_chrom, window_start, window_end))

            # plot the segmentation data
            for index, seg in seg_df.iterrows():
                for i, value in enumerate(self.segment_values):
                    start, end, plot_value = seg[['start', 'end', value]].values
                    ax.plot([start, end], [plot_value, plot_value], color=seg_colors[i],
                            **plot_kwargs)

        # handle axis limits depending on whether a window was passed
        if window is None:
            ax.set_xticks([loc for _, loc in xticks])
            ax.set_xticklabels([label for label, _ in xticks])
            ax.tick_params(axis='x', length=0)

            end = self.chrom_df['length'].sum()
            ax.set_xlim([0, end])
            ax.set_xlabel('')
        else:
            ax.set_xlim([window_start, window_end])

        if ylim is not None:
            ax.set_ylim(ylim)

        return ax

    def visualize_gene_profile(self, ax, gene, bed_df, padding=1e7, **kwargs):

        '''A simple wrapper function that identifies the start and stop of a gene
        from a bed file, then calls visualize_profile to visualize
        data around that gene.

        ax: axis object
            Axis object on which to plot the data
        gene: str
            name of gene to find in bed file
        bed_df: pandas dataframe
            Required columns are gene, chrom, start, and end
        padding: float
            region to expand around gene interval (default 1e6)
        kwargs: kwargs to be passed to visualize_coverage

        returns:
        -------
        ax: matplotlib axis object containing the figure
        '''

        required_cols = ['gene', 'chrom', 'start', 'end']
        missing_cols = set(required_cols) - set(bed_df.columns)
        if missing_cols:
            raise ValueError('Columns missing in bed file: {}'.format(missing_cols))

        gene_count = bed_df['gene'].to_list().count(gene)
        if gene_count == 0:
            raise ValueError('Gene not found in bed file')
        elif gene_count > 1:
            raise ValueError('Gene found in multiple entries in bed file')

        gene_index_bed_df = bed_df.set_index('gene')
        chrom, start, end = gene_index_bed_df.loc[gene, ['chrom','start', 'end']].values

        ax = self.visualize_profile(ax=ax, window=[chrom, start - padding, end + padding],
                                        **kwargs)

        ax.axvspan(start, end, color='grey', alpha=0.25)
        return ax
