import matplotlib.pyplot as plt
import pandas as pd
from sklearn.manifold import MDS


def get_ordination(sample_size):
    dist_matrix = pd.read_csv('KCL-Content/metadata/pw_dist_' + sample_size + '.csv', sep='\t')
    mds = MDS(dissimilarity='precomputed', random_state=0)
    ordination = mds.fit_transform(dist_matrix)
    pd.DataFrame(ordination).to_csv('KCL-Content/metadata/ordination_' + sample_size + '.csv', sep='\t', index=False)


def plot_ordination(var, sample_size):
    ordination = pd.read_csv('KCL-Content/metadata/ordination_' + sample_size + '.csv', sep='\t')
    pw_meta = pd.read_csv('KCL-Content/metadata/meta_' + sample_size + '.csv', sep='\t')

    fig, ax = plt.subplots()
    scatter = plt.scatter(ordination['0'], ordination['1'], s=2, c=pw_meta[var].astype('category').cat.codes)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.legend(handles=scatter.legend_elements()[0], 
            labels=list(set(pw_meta[var])),
            title=var)
    plt.savefig('KCL-Content/figures/MDS/' + sample_size + '_' + var + '.png', bbox_inches='tight') # tight = resize to labels

def plot_ordination_sizes(var, sample_size):
    ordination = pd.read_csv('KCL-Content/metadata/ordination_' + sample_size + '.csv', sep='\t')
    pw_meta = pd.read_csv('KCL-Content/metadata/meta_' + sample_size + '.csv', sep='\t')
    pw_meta['peptide'] = pw_meta['pool_id'].replace('A','', regex=True).replace('B','', regex=True).replace('Neg-','CTR', regex=True).replace('DMSO-','CTR', regex=True).replace('CEF-','CTR', regex=True)
    print(list(set(pw_meta['peptide'])))

    fig, ax = plt.subplots()
    scatter = plt.scatter(ordination['0'], ordination['1'], s=pw_meta['frequencyCount (%)']*500, c=pw_meta[var].astype('category').cat.codes)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.legend(handles=scatter.legend_elements()[0], 
            labels=list(set(pw_meta[var])),
            title=var)
    plt.savefig('KCL-Content/figures/MDS/' + sample_size + '_' + var + '_S.svg', bbox_inches='tight') # tight = resize to labels
    plt.savefig('KCL-Content/figures/MDS/' + sample_size + '_' + var + '_S.png', bbox_inches='tight') # tight = resize to labels


var='peptide'; sample_size='0.001'

# get_ordination(sample_size)
# plot_ordination(var, sample_size)
plot_ordination_sizes(var, sample_size)
