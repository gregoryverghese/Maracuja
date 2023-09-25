import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import circlify as circ
import seaborn as sns
from matplotlib.patches import Circle
from initialise_db import Patient, Sample, session, initialise_db
import matplotlib as mpl
from  matplotlib.colors import LinearSegmentedColormap

# pd.set_option('display.max_colwidth', None)
# pd.options.mode.chained_assignment = None # to silence warning about df self-edit

# initialise_db()


# -------- PLOTS ------------------------------------------------- 

def plot_bubbles(sample):   
    print(sample.id)
    df = sample.get_df()
    df2=df.sort_values('count', ascending = False).head(100)

    data = []
    for i in df2.index:
        data.append({'id': str(i), 'datum': df2['proportion'][i]})
    circles = circ.circlify(data)

    fig, ax = plt.subplots()
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    
    n = len(circles)
    for i in range(n):
        circle = circles[i]
        x = Circle((circle.x, circle.y), circle.r)
        ax.add_patch(x)
        # If clonotype expanded, label it:
        if circle.ex['datum']>0.05:
            plt.text(circle.x, circle.y, str(n - i))
    
    file_path = '../KCL-Content/figures/Sample_bubbles/'
    file_name = sample.id + '.svg'
    plt.savefig(file_path+file_name, bbox_inches='tight')
    plt.close()

# for sample in session.query(Sample).all():
#     plot_bubbles(sample)

# def expansion_heatmap():
metadata = pd.read_csv('../KCL-Content/metadata/metadata.csv', sep='\t')
expansion_index = metadata.pivot(index='patient-protocol',
                    columns='pool',
                    values='dominant_pro')
column_to_move = expansion_index.pop('Pool-10')
expansion_index.insert(9, 'Pool-10', column_to_move)
expansion_index = expansion_index.drop('Pool-T', axis=1).drop('KCL710-T', axis=0).drop('KCL725-T', axis=0).drop('KCL763-T', axis=0)
# print(expansion_index)
# sns.heatmap(expansion_index)
# plt.show()
# plt.savefig('../KCL-Content/figures/expansion_index/expansion_index-no-T.png', bbox_inches='tight')


# NO T (arch)
def one_plot_bubbles():
    plt.figure(figsize = (4,4))
    for sample in session.query(Sample).all():
        print(sample.id)
        df = sample.get_df()
        pool = sample.pool_id.replace('-A', '').replace('-B', '').replace('A', '').replace('B', '')
        patient_protocol = str(sample.patient_id) + '-' + str(sample.protocol)

        x_order = ['CEF', 'DMSO', 'Neg', 'Pool-1', 'Pool-2', 'Pool-5', 'Pool-7', 'Pool-8', 'Pool-9', 'Pool-10']
        y_order = ['KCL710-A', 'KCL710-B', 'KCL717-A', 'KCL725-A', 'KCL737-A', 'KCL750-A', 'KCL763-A', 'KCL763-B']

        if pool in x_order:
            x_row = x_order.index(pool)
            y_col = y_order.index(patient_protocol)

            data = []
            df2=df.sort_values('count', ascending = False).head(100)
            for i in df2.index:
                data.append({'id': str(i), 'datum': df2['count'][i]})
            circles = circ.circlify(data)

            ax = plt.subplot2grid((8,10), (y_col,x_row))
            ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
            ax.axis('off')
            
            sample_size = len(circles)
            for i in range(sample_size):
                circle = circles[i]
                x = Circle((circle.x, circle.y), circle.r)
                ax.add_patch(x)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('../KCL-Content/figures/expansion_index/bubble_overlay-non-T.svg', bbox_inches='tight')


# just A or B
def one_plot_bubbles():
    for sample in session.query(Sample).all():
        print(sample.id)
        df = sample.get_df()
        pool = sample.pool_id.replace('-A', '').replace('-B', '').replace('A', '').replace('B', '')
        patient_protocol = str(sample.patient_id) + '-' + str(sample.protocol)

        # expansion_index = df['count'].max()/63338
        expansion_index = df['proportion'].max()
        tcell_total = sum(df['count'])/138891

        x_order = ['CEF', 'DMSO', 'Neg', 'Pool-1', 'Pool-2', 'Pool-5', 'Pool-7', 'Pool-8', 'Pool-9', 'Pool-10']
        y_order = ['KCL710-A', 'KCL717-A', 'KCL725-A', 'KCL737-A', 'KCL750-A', 'KCL763-A', 'KCL710-B', 'KCL763-B']

        if 'T' not in sample.protocol:
            x_row = x_order.index(pool)
            y_col = y_order.index(patient_protocol)

            data = []
            df2=df.sort_values('count', ascending = False).head(100)
            for i in df2.index:
                data.append({'id': str(i), 'datum': df2['count'][i]})
            circles = circ.circlify(data)

            ax = plt.subplot2grid((8,10), (y_col,x_row))
            ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
            ax.axis('off')
            
            # colour = (tcell_total, 1 - tcell_total, 0)
            colour = (expansion_index, 1 - expansion_index, 0)

            sample_size = len(circles)
            for i in range(sample_size):
                circle = circles[i]
                x = Circle((circle.x, circle.y), circle.r, color=colour, linewidth=0)
                ax.add_patch(x)

    ax = plt.subplot2grid((9,10), (0,0))
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')

    for i in range(len(y_order)):
        plt.text(-5, -2.35*i, y_order[i], fontsize = 10)

    for i in range(len(x_order)):
        plt.text(2.45*i-0.5, -19.5, x_order[i], fontsize = 10, rotation=90)

    plt.text(9, -21, 'Treatment', fontsize = 10)
    plt.text(-6, -8, 'Patient', fontsize = 10, rotation=90)
    ax.axis('off')

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('../KCL-Content/figures/expansion_index/bubble_overlay_EI.svg', bbox_inches='tight')


one_plot_bubbles()

# df = session.query(Sample).first().get_df()
# print(df)
# cmap=LinearSegmentedColormap.from_list('rg',[(1,0,0), (0,1,0)], N=256)
# norm = mpl.colors.Normalize(vmin=5, vmax=10)
# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#             cax=ax, orientation='horizontal', label='Some Units')