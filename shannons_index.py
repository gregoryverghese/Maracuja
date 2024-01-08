import pandas as pd
import matplotlib.pyplot as plt
from database.initialise_db import Sample, session, initialise_db
import os
import math
import statistics
import scipy
import numpy as np
from sklearn.linear_model import LinearRegression
import sys

path = sys.argv[1] # '../KCL-Clean-data/'

initialise_db('../KCL-Clean-data')
samples = session.query(Sample).all()

metadata = pd.read_csv( path + 'metadata.csv', sep='\t')
shannons_stats = pd.DataFrame(index=metadata['Sample'])

### --- Add Shannons, Simpsons, and Gutais to metadata --- ###
if True:
    for sample in samples:
        df = sample.get_df()
        print(sample.id)

        # Shannon's = - Σ (proportion * ln proportion)
        def lnxpi(x):
            return x * math.log(x)
        pi = df['proportion']
        pilnpi = list(map(lnxpi, pi))
        shannons = round(- sum(pilnpi), 3)

        # Simpson's = 1 / Σ proportion **2
        def sq(x):
            return x**2
        pi = df['proportion']
        pi2 = list(map(sq, pi))
        simpsons = round(1/sum(pi2), 3)

        # Gutai's = proportion of dominant
        gutais = round(df['proportion'][0], 3)

        metadata.at[metadata[metadata['Sample'] == sample.id].index.values[0], 'shannons'] = shannons
        metadata.at[metadata[metadata['Sample'] == sample.id].index.values[0], 'simpsons'] = simpsons
        metadata.at[metadata[metadata['Sample'] == sample.id].index.values[0], 'gutais'] = gutais
        metadata.at[metadata[metadata['Sample'] == sample.id].index.values[0], 'N'] = sum(df['count'])
        
        shannons_stats.at[sample.id, 'shannons'] = shannons
        shannons_stats.at[sample.id, 'simpsons'] = simpsons
        shannons_stats.at[sample.id, 'gutais'] = gutais
        shannons_stats.at[sample.id, 'N'] = sum(df['count'])

    metadata.to_csv( path + 'metadata.csv', sep='\t', index=False)
    shannons_stats.to_csv( path + 'shannons.csv', sep='\t')


### Define function to statistically verify difference in Shannon's Index ###
if True:
    def getSS(x):
        return x * math.log(x) **2
    def getSQ(x):
        return x * math.log(x)

    def getStats(counts, log=False):
        pi = [ x / sum(counts) for x in counts ]
        N = sum(counts)
        richness = len(counts)

        SS = sum(list(map(getSS, pi)))
        SQ = sum(list(map(getSQ, pi))) ** 2
        H = abs(sum(list(map(getSQ, pi))))
        S2H = (SS-SQ)/N + (richness -1)/(2*N**2)

        if log:
            print('SS: ', SS)
            print('SQ: ', SQ)
            print('H: ', H)
            print('S2H: ', S2H)
            print("")
        return N, H, S2H

    # Compare 
    def compare_shannons(test, CTR):
        ctr_df = CTR.get_df()
        test_df = test.get_df()

        ctr_N, ctr_H, ctr_S2H = getStats(ctr_df['count'])
        test_N, test_H, test_S2H = getStats(test_df['count'])

        t = abs(test_H - ctr_H)/(test_S2H + ctr_S2H)**0.5
        df = round((test_S2H + ctr_S2H) **2 / ((test_S2H **2 /test_N) + (ctr_S2H **2 /ctr_N)), 0)
        p = scipy.stats.t.sf(abs(t), test_N + ctr_N -1)*2 
        return t, p

### Convert shannons.csv to treatment versus control ###
if True:
    shannons_stats = pd.read_csv( path + 'shannons.csv', sep='\t')

    treatment_samples = [ sample for sample in samples if 'Pool' in sample.id and 'T' not in sample.id ] 
    shannons_compare = pd.DataFrame(index=treatment_samples, columns=['shannons', 'DMSO', 'DMSO pval'])

    for treatment_sample in treatment_samples:
        print(treatment_sample.id)
        shannons_stat = shannons_stats['shannons'][shannons_stats['Sample']==treatment_sample.id].values[0]
        shannons_compare.at[treatment_sample, 'shannons'] = shannons_stat

        sample_name = treatment_sample.id.split('_')[0]

        for CTR_type in ['DMSO', 'Neg', 'CEF']:
            ctr_sample_name = sample_name + '_' +  CTR_type
            print(CTR_type)

            try:
                ctr_sample = [ sample for sample in samples if sample.id == ctr_sample_name ][0]

                shannons_stat = shannons_stats['shannons'][shannons_stats['Sample']==ctr_sample.id].values[0]
                print(shannons_stat)
                shannons_compare.at[treatment_sample, CTR_type] = shannons_stat
                t, p = compare_shannons(treatment_sample, ctr_sample)
                print(p)
                shannons_compare.at[treatment_sample, CTR_type + ' pval'] = p
            except IndexError:
                pass
            except ValueError:
                print("value error")
                
    print(shannons_compare)
    shannons_compare.to_csv( path + 'shannons_compare.csv', sep='\t')

def lnxpi(x):
    return x * math.log(x)
# Shannon's = - Σ (proportion * ln proportion)
def get_shannons(df):
    pi = df['proportion']
    pilnpi = list(map(lnxpi, pi))
    shannons = round(- sum(pilnpi), 3)
    return shannons

### Add to METADATA - treatment versus control ###
if True:
    metadata = pd.read_csv( path + 'metadata.csv', sep='\t', index_col=0)

    treatment_samples = [ sample for sample in samples if 'Pool' in sample.id and 'T' not in sample.id ] 

    for treatment_sample in treatment_samples:
        df = treatment_sample.get_df()
        print(treatment_sample.id)

        
        shannon = get_shannons(df)

        metadata.at[metadata[metadata.index == treatment_sample.id].index.values[0], 'shannons'] = shannon
        metadata.at[metadata[metadata.index == treatment_sample.id].index.values[0], 'N'] = sum(df['count'])
    
        sample_name = treatment_sample.id.split('_')[0]

        for CTR_type in ['DMSO', 'Neg', 'CEF']:
            ctr_sample_name = sample_name + '_' +  CTR_type

            try:
                ctr_sample = [ sample for sample in samples if sample.id == ctr_sample_name ][0]
                ctr_shannon = get_shannons(ctr_sample.get_df())
                metadata.at[treatment_sample.id, CTR_type + 'shannon'] = ctr_shannon
                t, p = compare_shannons(treatment_sample, ctr_sample)
                print(CTR_type)
                metadata.at[treatment_sample.id, CTR_type + 'pval'] = round(p, 3)
                metadata.at[treatment_sample.id, CTR_type + 'pval-bc55'] = round(p*55, 3)
            except IndexError:
                pass
            except ValueError:
                print("value error")
                
    print(metadata)
    metadata.to_csv( path + 'metadata.csv', sep='\t')

# # linear regression
if True:
    X = 'N'; Y = 'gutais'

    x = np.array(shannons_stats[X]).reshape((-1, 1))
    y = np.array(shannons_stats[Y])
    model = LinearRegression().fit(x, y)

    stat, p = scipy.stats.pearsonr(shannons_stats[X], shannons_stats[Y])
    p = round(p, 4)
    r_sq = round(model.score(x, y), 3)

    # plot
    plt.scatter(shannons_stats[X], shannons_stats[Y])

    x = np.linspace(0, 80000)
    y = model.coef_ * x + model.intercept_
    plt.plot(x, y)
    plt.xlabel(X); plt.ylabel(Y)
    plt.title(f"p: {p}, r sq = {r_sq}")
    plt.show()


