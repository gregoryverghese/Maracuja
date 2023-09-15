import pandas as pd 
import glob

pd.options.mode.chained_assignment = None # to silence warning about df self-edit

def format_1(data_path):
    df = pd.read_csv(data_path)

    for i in range(len(df[df.columns[0]])):
        df[df.columns[0]][i]=df[df.columns[0]][i].replace(df[df.columns[5]][i] + "_clonotype",'')

    for i in range(len(df[df.columns[5]])):
        df[df.columns[5]][i]=df[df.columns[5]][i].replace( "BIOKEY_",'')

    df.to_csv('processed_data/data1-BIOKEY.csv', sep=',')
    return(df)


def format_2(data_path):
    metadata = pd.read_csv(data_path + 'metadata.txt', delimiter='\t')
    df = pd.DataFrame()
    files = glob.glob(data_path + "/*.tsv")
    for file in files:
        file_name=file.replace(data_path,'').replace(".tsv",'')
        patient = metadata['patient_id'][metadata['Sample']==file_name].values[0]
        type = metadata['type'][metadata['Sample']==file_name].values[0]
        pool_id = metadata['pool_id'][metadata['Sample']==file_name].values[0]
        genomics = metadata['Genomics'][metadata['Sample']==file_name].values[0]
        hormad1 = metadata['HORMAD1'][metadata['Sample']==file_name].values[0]

        data = pd.read_csv(file, delimiter='\t')
        data = data.join(pd.DataFrame([patient] * len(data), columns=["patient"]))
        data = data.join(pd.DataFrame([type] * len(data), columns=["type"]))
        data = data.join(pd.DataFrame([pool_id] * len(data), columns=["pool_id"]))
        data = data.join(pd.DataFrame([genomics] * len(data), columns=["genomics"]))
        data = data.join(pd.DataFrame([hormad1] * len(data), columns=["hormad1"]))
        data = data.join(pd.DataFrame([file_name] * len(data), columns=["Sample"]))
        # df = df._append(data)
        df = pd.concat([df, data], ignore_index=True)

    df=df.rename(columns={"count (templates/reads)": "frequency"})
    df['clonotype_id'] = df.loc[:, 'nucleotide']
    df=df.set_axis(range(len(df)), axis=0)

    # df.to_csv('processed_data/data2-KCL.csv', sep=',')
    return df, metadata



# format_1('raw_data/1882-BIOKEY_clonotypes_combined_cohort2.csv')
# df2, metadata = format_2('raw_data/pbmc-v2/')
# print(df2)



# CHECK they are the same
# Has correct vers?
# print(len(df1['frequency']))
# print(len(df1['clonotype_id']))

