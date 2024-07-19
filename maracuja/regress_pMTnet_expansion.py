import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from database.initialise_db import Sample, session, initialise_db
import os
import math
import statistics
import scipy

# Get maximum predictions for sample for all peptides and HLA types

initialise_db('../KCL-Clean-data')
samples = session.query(Sample).all()

for sample in samples[0:1]:
    print(sample)
    path = 'pMTnet/sample_outputs'
    files = os.listdir(path)
    filess = [f.split('__')[0] for f in files]
    print(list(set(filess)))

    # files = [file for file in files if sample.id in file]
    # preds = pd.read_csv(path + '/' + files[0])[['CDR3', 'Rank']]
    # for file in files:
    #     next_pred = pd.read_csv(path + '/' + file)[['CDR3', 'Rank']].rename(columns = {'Rank':file})
    #     preds = pd.merge(preds, next_pred, on='CDR3').drop_duplicates()
    # max_preds = preds.min('columns')
    # total_preds = preds[['CDR3']]
    # total_preds['max_prediction'] = max_preds
    # preds = total_preds[['CDR3', 'max_prediction']].rename(columns = {'CDR3':'aminoacid'})


    # # Counts DF
    # counts = pd.read_csv('../KCL-Raw-data/' + sample + '.tsv', sep='\t')[['aminoAcid', 'count (templates/reads)']].rename(columns = {'count (templates/reads)':'count', 'aminoAcid':'aminoacid'})
    # counts = counts[counts['aminoacid'].notna()]
    # counts=counts.sort_values('count', ascending = False).reset_index().drop('index', axis=1)


    # # Merge predictions, ordination, and counts
    # # merged_df = ordination.merge(preds, on ='aminoacid',how='left').drop_duplicates()
    # # merged_df = merged_df.merge(counts, on ='aminoacid',how='left').drop_duplicates()
    # merged_df = counts.merge(preds, on ='aminoacid',how='left').drop_duplicates()

    # # Standardise between 0 and 1, create colour arrays
    # merged_df['max_prediction']=merged_df['max_prediction'].add(+9.9999999999989e-05).div(1.0001)
    # rank_col=[]
    # for rank in merged_df['max_prediction']:
    #     rank_col.append((1-rank, rank, 0))

    # count_col=[]
    # for count in merged_df['count']:
    #     proportion = count/max(merged_df['count'])
    #     count_col.append((proportion, 1 - proportion, 0))
    # # print(count_col)

    # # Sort for plot colour
    # merged_df = merged_df.sort_values(by='max_prediction', ascending=False)
    # merged_df = merged_df.sort_values(by='count', ascending=False)


    # # fig, ax = plt.subplots()
    # # # scatter = plt.scatter(merged_df['0'], merged_df['1'], s=1, c=count_col) 
    # # scatter = plt.scatter(merged_df['0'], merged_df['1'], s=3, c=rank_col) 
    # # ax.axis('off')
    # # plt.show()
    # # plt.savefig('../KCL-Content/figures/GIF/' + file + '.svg', format="svg")


    # # #### LINEAR REGRESSION
    # max_prediction = np.array(merged_df['max_prediction']).reshape((-1, 1))
    # counts = np.array(merged_df['count'])


    # model = LinearRegression()
    # model.fit(max_prediction, counts)
    # model = LinearRegression().fit(max_prediction, counts)
    # plt.scatter(max_prediction, counts)
    # plt.savefig('../KCL-Content/figures/pMTnet_count_regression/' + sample + '.png', format="png")

    # r_sq = "%.5f" % model.score(max_prediction, counts)
    # print(f"R^2 = {r_sq}")

    # max_prediction = np.array(merged_df['max_prediction'])
    # result = "%.5f" % np.corrcoef(max_prediction, counts)[0, 1]
    # print(f"p val = {r_sq}")

