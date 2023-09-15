

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import circlify as circ
from preprocess import format_1, format_2
from matplotlib.patches import Circle
import os
import re
pd.set_option('display.max_colwidth', None)
pd.options.mode.chained_assignment = None # to silence warning about df self-edit

BIOKEY = format_1('raw_data/1882-BIOKEY_clonotypes_combined_cohort2.csv')
KCL_bulk, metadata = format_2('raw_data/pbmc-v2/')

def get_expanded(df):
    total_tcells = sum(df['frequency'])
    proportion = df['frequency']/total_tcells
    expanded_nt = df['nucleotide'][proportion>0.05]
    expanded_aa = df['aminoAcid'][proportion>0.05]
    return expanded_aa, expanded_nt, proportion[proportion>0.05]


# T cell and clonotype total
def get_totals(df):
    tcell_total = sum(df['frequency'])
    unique_clonotypes = (list(set(df['clonotype_id'])))
    clonotype_total = len(unique_clonotypes)
    return tcell_total, clonotype_total


# -------- PLOTS ------------------------------------------------- 
# T cell and clonotype total per sample, respectively
def plot_tcell_clonotype_totals(df,var, file_path):
    plt.subplots(figsize=(20,20)) 
    factors = (list(set(df[var])))
    for factor in factors:
        tcell_total = sum(df['frequency'][df[var]==factor])
        clonotype_total = len(list(set(df['clonotype_id'][df[var]==factor])))
        plt.bar(factor, tcell_total, color='r')
        plt.bar(factor, clonotype_total, color='b')
    plt.xticks(rotation=90, fontweight='bold', fontsize='15')
    plt.yticks(fontweight='bold', fontsize='15')
    plt.xlabel('Sample', fontweight='bold', fontsize='15', horizontalalignment='center')
    plt.ylabel('T cell and clonotype count', fontweight='bold', fontsize='15', horizontalalignment='center')
    file_name = var + '.png'
    plt.savefig(file_path+file_name)
    plt.close()

# ------------------------------------------------------------------
def bar_plot_clonotype_proportions(df, var, file_path):
    factors = (list(set(df[var])))
    for factor in factors:
        x = df['frequency'][df[var]==factor].sort_values(ascending=False)
        x = x.cumsum().sort_values(ascending=False) / sum(x)
        for val in x:
            plt.bar(factor, val)
    plt.xticks(rotation=90, fontweight='bold', fontsize='8')
    plt.yticks(fontweight='bold', fontsize='8')
    plt.xlabel('Sample', fontweight='bold', fontsize='8', horizontalalignment='center')
    plt.ylabel('Clonotype proportion', fontweight='bold', fontsize='8', horizontalalignment='center')
    plt.subplots_adjust(bottom=0.2)
    file_name = 'barchart.png'
    plt.savefig(file_path+file_name)
    plt.close()

# Bubble chart plot of clonotype proportions by sample
def plot_bubbles(df):    
    file_path = 'KCL-Content/figures/Sample_bubbles/'
    sample = list(set(df['Sample']))[0]
    data = []
    df2=df.sort_values('frequency', ascending = False).head(100)
    expanded_aa, expanded_nt, proportion = get_expanded(df2)


    for i in df2['nucleotide'].index:
        data.append({'id': str(df2['nucleotide'][i]), 'datum': df2['frequency'][i]})
    circles = circ.circlify(data)

    fig, ax = plt.subplots()
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    
    sample_size = len(circles)
    for i in range(sample_size):
        circle = circles[i]
        x = Circle((circle.x, circle.y), circle.r)
        ax.add_patch(x)
        # If clonotype expanded, label it:
        if expanded_nt.str.contains(circle.ex['id']).any():
            plt.text(circle.x, circle.y, str(sample_size - i))
    
    file_name = sample + '.svg'
    plt.savefig(file_path+file_name, bbox_inches='tight')
    plt.close()
    return ",".join(expanded_aa.apply(str).values)


def do_all_bubbles():
    metadata['expanded'] = ['']*42
    samples = (list(set(KCL_bulk['Sample'])))
    for sample in samples:
        print(sample)
        df=KCL_bulk[KCL_bulk['Sample']==sample]
        x = plot_bubbles(df)
        metadata.at[metadata[metadata['Sample']==sample].index.values[0], 'expanded'] = x
    metadata.to_csv('KCL-Content/metadata/metadata.csv', sep='\t')


do_all_bubbles()