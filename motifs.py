import pandas as pd
import numpy as np
from palmotif import compute_pal_motif, svg_logo
from statistics import mode
import skunk
import os
import matplotlib.pyplot as plt

# POOL proportional
def aa_reads_by_pool_arrays(samples, folder):
    pools = list(set([sample.pool for sample in samples]))
    for pool in pools:
        print(pool)
        Pool_all = [sample for sample in samples if sample.pool == pool ]

        CDR3_tot = []
        for sample in Pool_all:
            df = sample.get_df()

            CDR3s = df[df['b_aminoacid'].notna()]
            CDR3s = df[~df['b_aminoacid'].str.contains('\*')]

            tot_count = sum(CDR3s['count'])

            for index, row in CDR3s.iterrows():
                aa = row['b_aminoacid']
                norm_count = round(row['count']/tot_count*10000)
                for i in range(norm_count):
                    CDR3_tot.append(aa)

        np.save(folder + '/' + pool, CDR3_tot)

# ALL proportional
def aa_reads_all_array(samples, folder):
    CDR3_tot = []
    for sample in samples:
        df = sample.get_df()

        CDR3s = df[df['b_aminoacid'].notna()]
        CDR3s = df[~df['b_aminoacid'].str.contains('\*')]

        tot_count = sum(CDR3s['count'])

        for index, row in CDR3s.iterrows():
            aa = row['b_aminoacid']
            norm_count = round(row['count']/tot_count*10000)
            for i in range(norm_count):
                CDR3_tot.append(aa)
    np.save(folder + '/all', CDR3_tot)

# SAMPLE proportional
def aa_reads_by_sample_arrays(samples, folder):
    for sample in samples:
        CDR3_tot = []
        df = sample.get_df()

        CDR3s = df[df['b_aminoacid'].notna()]
        CDR3s = df[~df['b_aminoacid'].str.contains('\*')]

        tot_count = sum(CDR3s['count'])

        for index, row in CDR3s.iterrows():
            aa = row['b_aminoacid']
            norm_count = round(row['count']/tot_count*10000)
            for i in range(norm_count):
                CDR3_tot.append(aa)
        
        np.save(folder + '/' + sample.id, CDR3_tot)

# # CTRs proportional
def aa_reads_CTR_array(samples, folder):
    CTRs = [sample for sample in samples if 'Pool' not in sample.pool]
    CDR3_tot = []
    for sample in CTRs:
        df = sample.get_df()

        CDR3s = df[df['b_aminoacid'].notna()]
        CDR3s = df[~df['b_aminoacid'].str.contains('\*')]

        tot_count = sum(CDR3s['count'])

        for index, row in CDR3s.iterrows():
            aa = row['b_aminoacid']
            norm_count = round(row['count']/tot_count*10000)
            for i in range(norm_count):
                CDR3_tot.append(aa)
    np.save(folder + '/CTRs', CDR3_tot)

# # SAMPLE each
def aa_by_sample_arrays(samples, folder):
    for sample in samples:
        CDR3_tot = []
        df = sample.get_df()

        CDR3s = df[df['b_aminoacid'].notna()]
        CDR3s = CDR3s[~CDR3s['b_aminoacid'].str.contains('\*')]
        CDR3s = CDR3s.head(20)
        CDR3_tot = list(CDR3s['b_aminoacid'])

        np.save(folder + '/' + sample.id, CDR3_tot)


# # POOL each
def aa_by_pool_arrays(samples, folder):
    pools = list(set([sample.pool for sample in samples]))
    for pool in pools:
        Pool_all = [sample for sample in samples if sample.pool == pool ]

        CDR3_tot = []
        for sample in Pool_all:
            df = sample.get_df()

            CDR3s = df[df['b_aminoacid'].notna()]
            CDR3s = df[~df['b_aminoacid'].str.contains('\*')]
            CDR3s = CDR3s.head(20)
            CDR3_tot = list(CDR3s['b_aminoacid'])

        np.save(folder + '/' + pool, CDR3_tot)

def arrays_to_logos(input_dir, output_dir):
    files = os.listdir(input_dir)
    for file in files:
        file_name = os.path.splitext(file)[0]
        print(file_name)
        CDR3_tot = np.load(input_dir + '/' + file)

        centroid = mode(CDR3_tot)

        motif, stats = \
            compute_pal_motif(
                centroid = centroid,
                seqs = CDR3_tot,
            )

        svg = svg_logo(motif, return_str = True)
        svg_logo(motif, output_dir + '/' + file_name + '.svg')


def arrays_to_logos_refs(input_dir, output_dir):
    files = os.listdir(input_dir)
    for file in files:
        file_name = os.path.splitext(file)[0]
        print(file)
        CDR3_tot = np.load(input_dir + '/' + file)
        centroid = mode(CDR3_tot)
        refs = np.load(input_dir + '/CTRs.npy')

        motif, stats = \
            compute_pal_motif(
                centroid = centroid,
                seqs = CDR3_tot,
                refs = refs
            )

        svg = svg_logo(motif, return_str = True)
        svg_logo(motif, output_dir + '/' + file_name + '.svg')

def plot_motifs(logos_folder, x_order, y_order, plot_folder):
    files = os.listdir(logos_folder)

    for file in [file for file in files if 'KCL' in file]:

        file_name = os.path.splitext(file)[0]

        patient_protocol = file_name.split('_')[0]
        pool = file_name.split('_')[1]

        x_row = x_order.index(pool)
        y_col = y_order.index(patient_protocol)

        ax = plt.subplot2grid((len(y_order) + 1, len(x_order)), (y_col,x_row))
        ax.set_xlim(0, 1); ax.set_ylim(0.3, 0.7); ax.set_aspect('equal')
        ax.axis('off')

        skunk.connect(ax, file_name)
    plt.subplots_adjust(wspace=0, hspace=0)

    for i in range(len(y_order)):
        plt.text(-4, 0.83*i-2.9, y_order[i], fontsize = 5)
    for i in range(len(x_order)):
        plt.text(i -2.5, -4.5, x_order[i], fontsize = 5, rotation=90)

    for file_name in x_order:

        file = file_name + '.svg'
        x_row = x_order.index(file_name)
        y_col = len(y_order) 

        ax = plt.subplot2grid((len(y_order) + 1, len(x_order)), (y_col,x_row))
        ax.set_xlim(0, 1); ax.set_ylim(0.3, 0.7); ax.set_aspect('equal')
        ax.axis('off')

        skunk.connect(ax, file_name)
    plt.subplots_adjust(wspace=0, hspace=0)


    dicts = {}
    keys = [os.path.splitext(file)[0] for file in files]
    values = [file for file in files]
    for i in range(len(keys)):
        dicts[keys[i]] = logos_folder + '/' + values[i]

    print(dicts)
    svg = skunk.insert(dicts)
    # svg = skunk.insert({'KCL737-A_Pool1': '../KCL-content/arrays/logos/KCL737-A_Pool1.svg'})
    

    with open(plot_folder + '/proportional_samples1.svg', 'w') as f:
        f.write(svg)

arrays_folder = 'arrays/proportional'
logos_folder = 'logos'
arrays_to_logos_refs(arrays_folder, logos_folder)


### MOTIF ANALYSIS
arrays_folder = '../KCL-content/arrays'
logos_folder = '../KCL-content/logos'
plot_folder = '../KCL-content/figures/motif_analysis'
aa_reads_CTR_array(samples, arrays_folder)
arrays_to_logos_refs(arrays_folder, logos_folder) # To run in docker container if Palmotif is not able to install on your local
plot_motifs('../KCL-content/refs_logos_each', x_order, y_order, plot_folder)


## proportional
motifs1.aa_reads_by_pool_arrays(samples, '../KCL-content/arrays/proportional')
create_arrays.aa_reads_by_sample_arrays(samples, '../KCL-content/arrays/proportional')
create_arrays.aa_reads_CTR_array(samples, '../KCL-content/arrays/proportional')

logos_to_fig.plot_motifs('../KCL-content/arrays/logos', x_order, y_order, '../KCL-content/figures/motif_analysis')
