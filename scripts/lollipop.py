import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable
import seaborn as sns
import requests
import json
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec


class LollipopPlot:

    '''Class for handling lolliplot plots for individual genes
    
    Attributes:
    -----------
    mutations: pandas dataframe
        Mutations in maf form. Relevant columns are as follows:
        
        Hugo_Symbol: gene name. An error will be raised if multiple
        genes are present

        Variant_Classification: Type of mutation (all will be shown)

        value_col: The height of the lollipop. This could be anything, from
        tumor fraction, to number of mutations

        Location: Genomic location of the alteration in the region/gene. 
        For genes, this should be amino acid #.

    gene: Name of gene (Hugo_Symbol) to plot.

    value_col: the column in the mutation_df that contains the 
        value that should be plotted

    length: int
        Length of the gene/region

    domains: list-like of list-like
        List containing domain information. Entries correspond
        to name, start, end, and color. Eg

        domains = [['P53 binding', 45, 64, 'blue'],
                   ['Promoter', 1, 16, 'yellow']] 
    '''

    def __init__(self, mutation_df, gene, value_col, length=None, domains=None):

        self.mutations = mutation_df

        required_cols = set(['Hugo_Symbol', 'Variant_Classification', value_col, 'Location'])
        missing_cols = required_cols - set(mutation_df.columns)
        if missing_cols:
            raise ValueError('Missing columns: {}'.format(missing_cols))

        self.value_col = value_col
        self.gene = gene
        self.length = length

        if domains is None:
            self.domains = []

    def load_gene_info(self, colors=None):

        '''Function that loads uses GenomeNexus and Pfam to fetch
        gene length and relevant domains. Sometimes fails due to
        what I presume are connection errors to GenomeNexus.

        Parameters:
        -----------
        colors: list of mpl.colors
            Define the color of domains (from left to right). Defaults
            to the current color cycle.

        Returns:
        --------
        pfamDomains:
            list of pfam domains that are associated with the gene
        '''

        request_url = 'https://www.genomenexus.org/ensembl/canonical-transcript/hgnc/{}?isoformOverrideSource=uniprot'.format(self.gene)
        result = requests.get(request_url)
        gene_data = json.loads(result.text)

        self.length = gene_data['proteinLength']

        pfamDomains = gene_data['pfamDomains']
        unique_domains = set([dom['pfamDomainId'] for dom in pfamDomains])

        if colors is None:
            default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            colors = [default_colors[i % len(default_colors)] for i in range(len(unique_domains))]

        if len(colors) < len(unique_domains):
            raise ValueError('Only {} colors defined for {} unique domains.'.format(len(colors)), len(unique_domains))

        # transform color to an iterator
        colors = iter(colors)

        # get domains, use a dict for caching
        domain_cache = {}

        sorted_domains = sorted(pfamDomains, key = lambda x: x['pfamDomainStart'])
        for pfam_dom in sorted_domains:

            # get the domain name
            pfam_id = pfam_dom['pfamDomainId']
            if pfam_id in domain_cache:
                name, col = domain_cache[pfam_id]
            else:
                pfam_request = 'https://www.genomenexus.org/pfam/domain/{}'.format(pfam_id)
                dom_data = json.loads(requests.get(pfam_request).text)
                name = dom_data['name']
                col = next(colors)

            # get start, end, and color
            start, end = pfam_dom['pfamDomainStart'], pfam_dom['pfamDomainEnd']

            # update domain dict with name and color
            pfam_dom['name'] = name

            # update cache
            domain_cache[pfam_id] = [name, col]

            self.domains.append([name, start, end, col])

        return pfamDomains

    def plot_lollipop(self, fig=None, spec=None, padding=None, figsize = (18,6), ylim = None,
                      height_ratios = (5,1), hspace=0, **scatter_kwargs):

        '''Plots the lollipop plot. Returns fig, (mut_ax, gene_ax)'''

        # get defaults
        if ylim is None:
            ylim = [0, None]

        if padding is None:
            padding = self.length // 50

        lollipop_data = self.mutations[self.mutations['Hugo_Symbol'] == self.gene].copy()

        if fig is None:
            fig = plt.figure(figsize=figsize)

        if spec is None:
            spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig,
                                     height_ratios=height_ratios,
                                     hspace=hspace)

        # otherwise, create gridspec in spec
        else:
            spec = gridspec.GridSpecFromSubplotSpec(ncols=1, nrows=2,
                                                    height_ratios=height_ratios,
                                                    hspace=hspace, subplot_spec=spec)

        mut_ax, gene_ax = fig.add_subplot(spec[0]), fig.add_subplot(spec[1])

        # plot the gene axis
        gene_ax.set_ylim([0,1])
        gene_ax.set_xlim([-padding, self.length + padding]) 
        gene_ax.add_patch(patches.Rectangle((0, 0.25), self.length, 0.5, color = 'grey', alpha = 0.5)) 
        gene_ax.set_xlabel('amino acid')

        # plot the domains, in the order 
        legend_patches = {}
        sorted_domains = sorted(self.domains, key = lambda x: x[1])
        for dom in sorted_domains:
            name, start, end, color = dom
            patch = patches.Rectangle((start, 0.15), end - start, 0.70, color = color, label = name)

            if name not in legend_patches:
                legend_patches[name] = patch

            gene_ax.add_patch(patch)

        # add gene legend, aligning title to the left
        gene_legend = gene_ax.legend(handles = legend_patches.values(), ncol = 1, 
                                     loc = 'lower left', bbox_to_anchor = (1, 0), title = 'domain')
        gene_legend.set_title('Domains')
        gene_legend._legend_box.align = "left"

        # clip the spine on the bottom axis
        gene_ax.spines['bottom'].set_bounds(0, self.length)

        # make sure the final tick mark is always present
        locs = list(gene_ax.get_xticks())
        locs = [int(loc) for loc in locs if (loc >= 0) and (loc <= self.length)] + [self.length]
        labels = [str(loc) for loc in locs]

        gene_ax.set_xticks(locs)
        gene_ax.set_xticklabels(labels)

        # remove the top, left, and right bounds
        gene_ax.get_yaxis().set_visible(False)
        sns.despine(ax = gene_ax, left = True)

        # plot mutation axis - vlines first
        mut_ax.vlines(x = lollipop_data['Location'], ymin = 0, ymax = lollipop_data[self.value_col], 
                      color='black')

        # handle defaults for scatter
        if 's' not in scatter_kwargs:
            scatter_kwargs['s'] = 120

        sns.scatterplot(data = lollipop_data, x = 'Location', y = self.value_col, hue = 'Variant_Classification', ax = mut_ax,
                        zorder = 2, **scatter_kwargs)

        # connect with the gene axis - kind of hacky. First we determine if any line is on a domain:
        lollipop_data['domain?'] = lollipop_data['Location'].apply(lambda loc: any([(loc > dom[1]) and (loc < dom[2]) for dom in self.domains]))
        ymin = 0.75 + 0.1*lollipop_data['domain?']

        # then we connect the mutation axis to the gene axis
        gene_ax.vlines(x = lollipop_data['Location'], ymin = ymin, ymax = 1, color = 'black')

        # alter mutation axis
        mut_ax.set_ylim(ylim)
        mut_ax.tick_params(axis = 'y')
        sns.despine(ax = mut_ax, bottom = True)

        # make mutation legend
        handles, labels = mut_ax.get_legend_handles_labels()

        mut_leg = mut_ax.legend(handles=handles, labels=labels, bbox_to_anchor = (1, 0.45), handletextpad = 0.7)
        mut_leg.set_title('Mutation type')
        for handle in mut_leg.legendHandles:
            handle._sizes = [scatter_kwargs['s']]

        mut_ax.set_xlim(gene_ax.get_xlim())

        return fig, (mut_ax, gene_ax)










