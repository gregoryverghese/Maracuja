import pandas as pd
import os
import sys

# Use this script to convert Raw-data in format stored in KCL-Clean-data

format = sys.argv[1]
path = sys.argv[2]

### --- Second lot of bulkDNAseq data (from immunoseq adaptive)--- ###
if int(format) == 2:
    # path = '../ARCHIVE/KCL-New-data'
    files = [ file for file in os.listdir(path) if 'KCL' in file ]

    for file in files:
        print(file)
        df = pd.read_csv(path + '/' + file, delimiter='\t')

        productive_templates = df['productive_templates'][0]
        df['productive_frequency'] = df['productive_frequency']*productive_templates
        df = df[~df['productive_frequency'].isna()]
        df = df.rename(columns={'rearrangement': 'nucleotide', 'amino_acid': 'aminoacid', 'productive_frequency': 'count'})
        df = df[['nucleotide', 'aminoacid', 'count']]
        df = df.sort_values('count', ascending=False)
        df['proportion'] = df['count']/sum(df['count'])

        df['count'] = df['count'].round(0).astype(int)
        df.to_csv('../KCL-Clean-data/' + file, sep='\t')


### --- First lot of bulkDNAseq data (from Greg) --- ###
if int(format) == 1:
    # path = '../KCL-Raw-data'
    files = [ file for file in os.listdir(path) if 'KCL' in file ]
    for file in files:
        print(file)
        df = pd.read_csv(path + '/' + file, delimiter='\t')
        df = df[['nucleotide', 'aminoAcid', 'count (templates/reads)']]
        df = df.rename(columns={'aminoAcid': 'aminoacid', 'count (templates/reads)': 'count'})
        df = df[df['aminoacid'].notna()].sort_values('count', ascending=False)
        
        df['proportion'] = df['count']/sum(df['count'])
        df.to_csv('../KCL-Clean-data/' + file, sep='\t', index = False)

