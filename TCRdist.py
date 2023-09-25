import pandas as pd
import pwseqdist as pw
from tcrdist.rep_funcs import _pws
import matplotlib.pyplot as plt
import numpy as np

# PARAMETERS:
threshold = 0.1 # threshold for clonotype expansion
closeness_threshold = 30 # threshold for nucleotide similarity


# Plot heatmap of TCR distances for AminoAcid in given dataframe
def get_pw_dist(df):
        dmats = _pws(df = df,
                metrics = { 'aminoAcid' : pw.metrics.nb_vector_tcrdist}, 
                weights= { 'aminoAcid' : 3}, 
                kargs={ 'aminoAcid' : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False}}, 
                cpu = 1, 
                store = True)

        pw_dist = dmats['tcrdist']
        # print(pw_dist)
        # print(df['aminoAcid'].to_list())
        return pw_dist

def plot_heatmaps(pw_dist, sample):
        n = len(pw_dist)
        fig, ax1 = plt.subplots()
        heatmap = ax1.pcolor(pw_dist, cmap='hot')

        plt.colorbar(heatmap)
        ax2 = ax1.twiny() # copy to define 2 axes
        ax2.axes.get_xaxis().set_visible(False) # hide top x labs
        ax1.set_xticks(np.linspace(0.5, n-0.5, num=n))
        ax1.set_xticklabels(dfa['aminoAcid'].to_list(), rotation=90) # label bottom x axis
        ax2.set_yticks(np.linspace(0.5, n-0.5, num=n))
        ax2.set_yticklabels(dfa['aminoAcid'].to_list()) # label left x axis

        plt.savefig('Analysis/KCL-Content/figures/Sample_expanded_heatmaps/' + sample + '.png', bbox_inches='tight') # tight = resize to labels
        plt.close()

def archive_pw_dist(sample_size):
        df = KCL_bulk[KCL_bulk['aminoAcid'].notna()] # remove nan
        df2 = df[df['frequencyCount (%)']>sample_size]
        # df2 = df.sample(n=500)
        dfa = df2.drop_duplicates(subset=["aminoAcid"], keep='first') #remove duplicates

        pw_dist = get_pw_dist(dfa)
        pd.DataFrame(pw_dist).to_csv('Analysis/pw_dist_' + str(sample_size) + '.csv', sep='\t', index=False)
        pd.DataFrame(dfa).to_csv('Analysis/meta_' + str(sample_size) + '.csv', sep='\t', index=False)


KCL_bulk = pd.read_csv("Analysis/processed_data/data2-KCL.csv", delimiter=',')

archive_pw_dist(0.001)

# # Plot per sample
# samples = (list(set(KCL_bulk['Sample'])))
# for sample in samples:
#         print(sample)
#         df = KCL_bulk[KCL_bulk['Sample']==sample]
#         df2 = df[df['aminoAcid'].notna()] # remove nan
#         dfa = df2.sort_values('frequency', ascending = False).head(10)
#         pw_dist = get_pw_dist(dfa)
#         plot_heatmaps(pw_dist, sample)



# Plot ALL



# plot_heatmaps(pw_dist, 'ALL_expanded')

# print(np.where((pw_dist > 0) & (pw_dist < 40)))
