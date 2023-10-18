import skunk
import numpy as np
import os
import matplotlib.pyplot as plt

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

    svg = skunk.insert(dicts)

    with open(plot_folder + '/all_samples.svg', 'w') as f:
        f.write(svg)
