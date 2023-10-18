import pandas as pd
import numpy as np
from palmotif import compute_pal_motif, svg_logo
from statistics import mode
import os


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


arrays_folder = '../KCL-content/arrays'
logos_folder = '../KCL-content/logos'
arrays_to_logos_refs(arrays_folder, logos_folder)