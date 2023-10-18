import pandas as pd
import numpy as np

# POOL
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

# ALL
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

# EACH
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

# # CTRs
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



# # SAMPLE EACH
def aa_by_sample_arrays(samples, folder):
    for sample in samples:
        CDR3_tot = []
        df = sample.get_df()

        CDR3s = df[df['b_aminoacid'].notna()]
        CDR3s = CDR3s[~CDR3s['b_aminoacid'].str.contains('\*')]
        CDR3s = CDR3s.head(20)
        CDR3_tot = list(CDR3s['b_aminoacid'])

        np.save(folder + '/' + sample.id, CDR3_tot)


# # POOL EACH
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