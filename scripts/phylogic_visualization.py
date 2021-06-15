import pandas as pd
import numpy as np
import seaborn as sns

class PhylogicTwoTPResult:
    
    def __init__(self, cluster_ccf_df, mutation_ccf_df):
        
        '''A class for visualizing two timepoint/biopsy Phylogic results'''
        
        self.cluster_ccf_df = cluster_ccf_df
        self.mutation_ccf_df = mutation_ccf_df
        self.biopsy_sites = ()
    
    def plot_phylogic(self, ax, color_palette = None, size_cutoff = 7, ccf_cutoff = 0.10):
        '''Plots phylogic result. Color palette is list of colors. If a cluster has fewer mutations than size_cutoff, it is ignored.
        If a cluster has ccf below ccf_cutoff for all biopsies, it is ignored.'''

        # calculate cluster sizes
        cluster_sizes = self.mutation_ccf_df.groupby('Cluster_Assignment').size()/2
        
        for i, cluster_id in enumerate(self.cluster_ccf_df['Cluster_ID'].drop_duplicates()):
            cluster_size = int(cluster_sizes[cluster_id])

            if cluster_size < size_cutoff:
                continue
    
            cluster_data = self.cluster_ccf_df[self.cluster_ccf_df['Cluster_ID'] == cluster_id]
        
            # dont show clusters with <= 0.10 ccf across all biopsies
            if all(cluster_data['postDP_ccf_mean'] < ccf_cutoff):
                continue
            
            if len(cluster_data) != 2:
                raise ValueError('Should be two timepoints')
            
            if color_palette is None:
                # plot mean ccf
                ax.plot([1,0], cluster_data['postDP_ccf_mean'], label = cluster_size)

                # plot confidence interval
                ax.fill_between(x = [1,0], y1 = cluster_data['postDP_ccf_CI_low'], 
                                y2 = cluster_data['postDP_ccf_CI_high'], alpha = 0.5)
            else:
                cluster_color = color_palette[i]

                # plot mean ccf
                ax.plot([1,0], cluster_data['postDP_ccf_mean'], label = cluster_size, color = cluster_color)

                # plot confidence interval
                ax.fill_between(x = [1,0], y1 = cluster_data['postDP_ccf_CI_low'], 
                                y2 = cluster_data['postDP_ccf_CI_high'], alpha = 0.5, color = cluster_color)

        ax.legend(bbox_to_anchor = (1,0.5), fontsize = 10, title = 'mutations')
        ax.set_xticks([0,1])
        ax.set_xlim([-0.08, None])
        bbx_label = str(self.biopsy_sites[0])
        ffpe_label = str(self.biopsy_sites[1])
        ax.set_xticklabels([ffpe_label, bbx_label], fontsize = 12)
        ax.tick_params(axis = 'x', which = 'both', length = 0)
        ax.set_ylabel('cancer cell fraction', fontsize = 10)

        sns.despine(bottom = True)
        
        return ax