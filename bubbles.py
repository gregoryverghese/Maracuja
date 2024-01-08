import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import circlify as circ
from matplotlib.patches import Circle
from database.initialise_db import Sample, session, initialise_db
import sys

path = sys.argv[1] # '../KCL-Clean-data'
outdir = sys.argv[2] # '../KCL-Content/figures/Sample_bubbles/'
initialise_db(path) 

# CREATE BUBBLES IMAGE FILES FOR ALL COUNT FILES IN DIRECTORY

def bubble(sample):
    print(sample)
    df = sample.get_df().head(100)

    data = []
    for i in df.index:
        data.append({'id': str(i), 'datum': df['count'][i]})
    circles = circ.circlify(data)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
    ax.axis('off')
    sample_size = len(circles)
    for i in range(sample_size):
        circle = circles[i]
        x = Circle((circle.x, circle.y), circle.r, linewidth=0)
        ax.add_patch(x)

samples = session.query(Sample).all()
for sample in samples:
    bubble(sample)
    plt.savefig(outdir + sample.id + '.svg', format='svg')
    plt.close()


### BUBBLE OVERLAY (no longer in use as we now have more, too many, samples)

def get_CTR(samples):
    CTRs = [sample for sample in samples if 'Pool' not in sample.pool]
    CTR = []
    for control in CTRs:
        df = control.get_df()
        aas = list(df['b_aminoacid'])
        CTR = CTR + aas
    return CTR

def bubble_overlay(samples, x_order, y_order, factor):

    for sample in samples:
        print(sample)
        df = sample.get_df()
        patient_protocol = str(sample.patient_id) + '-' + str(sample.protocol)

        x_row = x_order.index(sample.pool)
        y_col = y_order.index(patient_protocol)

        ax = plt.subplot2grid((len(y_order),len(x_order)), (y_col,x_row))
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
        ax.axis('off')

        colour = []
        if (factor == 'none'):
            expansion_index = df['proportion'].max()
            colour = [(0.118, 0.506, 0.690)]*len(df)

        if (factor == 'EI'):
            expansion_index = df['proportion'].max()
            EI = [expansion_index]*len(df)
            colour = [(x,1-x,0) for x in EI]

        if (factor == 'TCT'):
            tcell_total = sum(df['count'])/73866
            TCT = [tcell_total]*len(df)
            colour = [(x,1-x,0) for x in TCT]

        if (factor == 'in_control'):
            CTR = get_CTR(samples)
            for aa in df['b_aminoacid']:
                IC = int(aa in CTR)
                colour.append((0.118+IC*0.6, 0.506+IC*0.4, 0.690+IC*0.31)) 

        if (factor == 'pMTnet'):
            PT = df['pMTnet']
            colour = [(x,1-x,0) for x in PT]

        df['colour'] = colour

        data = []
        for i in df.index:
            data.append({'id': str(i), 'datum': df['count'][i], 'col': df['colour'][i]})
        circles = circ.circlify(data)

        sample_size = len(circles)
        for i in range(sample_size):
            circle = circles[i]
            colour = circle.ex['col']
            x = Circle((circle.x, circle.y), circle.r, color=colour, linewidth=0)
            ax.add_patch(x)

    ax = plt.subplot2grid((len(y_order),len(x_order)), (0,0))
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')

    for i in range(len(y_order)):
        plt.text(-6, -2.05*i, y_order[i], fontsize = 10)
    for i in range(len(x_order)):
        plt.text(2.48*i-0.5, -24, x_order[i], fontsize = 10, rotation=90)

    plt.text(9, -26.5, 'Treatment', fontsize = 10)
    plt.text(-7, -8, 'Patient', fontsize = 10, rotation=90)    
    ax.axis('off')

    plt.subplots_adjust(wspace=0, hspace=0)


# samples = [sample for sample in session.query(Sample).all() if sample.protocol != 'T' ]
# samples = [sample for sample in samples if '184+187' not in sample.id ]
# samples = [sample for sample in samples if 'AA' not in sample.id ]

# x_order = ['CEF', 'DMSO', 'Neg', 'Pool1', 'Pool2', 'Pool3', 'Pool5', 'Pool6', 'Pool7', 'Pool8', 'Pool9', 'Pool10']#, 'PoolT']
# y_order = ['KCL710-A', 'KCL717-A', 'KCL725-A', 'KCL737-A', 'KCL750-A', 'KCL763-A', 'KCL710-B', 'KCL763-B', 'KCL751-B', 'KCL1148-B', 'KCL1077-B']#, 'KCL710-T', 'KCL725-T', 'KCL763-T']
# bubble_overlay(samples, x_order, y_order, 'EI')
# plt.savefig('bubble_overlay_EI.svg', format='svg')

