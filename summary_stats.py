import numpy as np
import pandas as pd
import scipy.stats.distributions as dist
import os
from scipy.stats import fisher_exact
import sys

path = sys.argv[1] # '../KCL-Clean-data/'
treatment_files = [ file for file in os.listdir(path) if 'Pool' in file ]
df = pd.DataFrame()


# Fisher's Exact test (split by file)
for file in treatment_files[0:0]: 
    test = pd.read_csv(path + '/' + file, delimiter='\t')
    sample_name = file.split('_')[0]
    treatment = file.split('.')[0]

    print(treatment)
    df = pd.DataFrame()
    for i in range(5):
        clonotype = test['aminoacid'][i]

        test_count = test['count'][i]
        n_test = sum(test['count'])
        test_rest = (n_test - test_count)

        for CTR_type in ['DMSO', 'Neg', 'CEF']:
            CTR_file = path + '/' + sample_name + '_' + CTR_type + '.tsv'
            if os.path.isfile(CTR_file):            
                CTR = pd.read_csv(path + '/' + CTR_file, delimiter='\t')
            
                n_CTR = sum(CTR['count'])
                try:
                    CTR_count = CTR['count'][CTR['aminoacid']==clonotype].values[0]
                except IndexError:
                    CTR_count = 0
                CTR_rest = (n_CTR - CTR_count)

                table = pd.DataFrame({'test':[test_count, test_rest], 'control':[CTR_count, CTR_rest]}, index=pd.Index(['clonotype count', 'else count']))
                oddsr, p = fisher_exact(table=table.to_numpy(), alternative='two-sided')

                df.at[clonotype, 'count'] = round(test_count, 0)
                df.at[clonotype, 'proportion'] = round(test_count/n_test, 2)

                if (CTR_count/n_CTR) > (test_count/n_test):
                    df.at[clonotype, CTR_type] = '..'
                elif p<0.0005:
                    df.at[clonotype, CTR_type] = '***'
                elif p<0.005:
                    df.at[clonotype, CTR_type] = '**'
                elif p<0.05:
                    df.at[clonotype, CTR_type] = '*'
                else:
                    df.at[clonotype, CTR_type] = '.'

    df = df.reset_index().rename(columns={'index':'clonotype'})

    df.to_csv('../summary_stats/' + file, sep='\t', index=False)


# Control files
CTR_files = [ file for file in os.listdir(path) if 'Pool' not in file ]

for file in CTR_files: 
    df = pd.DataFrame()
    ctr = pd.read_csv(path + '/' + file, delimiter='\t')
    treatment = file.split('.')[0]

    for i in range(5):
        clonotype = ctr['aminoacid'][i]
        count = ctr['count'][i]
        total = sum(ctr['count'])

        df.at[clonotype, 'count'] = round(count, 0)
        df.at[clonotype, 'proportion'] = round(count/total, 2)

    df = df.reset_index().rename(columns={'index':'clonotype'})
    df.to_csv('../summary_stats/' + file, sep='\t', index=False)



# Chi-Squared test (not actually appropriate for samples with many zeros)
if False:
    test = pd.read_csv('../KCL-Clean-data/KCL763-B_Pool9.tsv', delimiter='\t')
    DMSO = pd.read_csv('../KCL-Clean-data/KCL763-B_DMSO.tsv', delimiter='\t')
    clonotype = test['aminoacid'][0]

    test_count = test['count'][0] # 21
    DMSO_count = DMSO['count'][DMSO['aminoacid']==clonotype].values[0] # 2

    test_proportion = test['proportion'][0] # 0.09292
    DMSO_proportion = DMSO['proportion'][DMSO['aminoacid']==clonotype].values[0] # 0.00029

    n_test = sum(test['count']) # 226
    n_DMSO = sum(DMSO['count']) # 6699
    mean_proportion = (test_count + DMSO_count)/(n_test + n_DMSO)

    # Standard error
    variance = mean_proportion * (1 - mean_proportion)
    standard_error = np.sqrt(variance * (1 / DMSO_count + 1 / DMSO_count))


    # Test statistic 
    best_estimate = (test_proportion - DMSO_proportion)
    hypothesized_estimate = 0
    test_stat = (best_estimate-hypothesized_estimate) / standard_error


    # P-value
    pvalue = 2*dist.norm.cdf(-np.abs(test_stat)) # Multiplied by two indicates a two tailed testing.
    pvalue = dist.norm.cdf(-np.abs(test_stat)) # Multiplied by one indicates a two tailed testing.
    print("P-value:", round(pvalue, 5))
