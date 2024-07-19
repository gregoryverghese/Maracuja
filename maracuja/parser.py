import pandas as pd
import os
from typing import Dict, Tuple


def is_immunoseq(file_path: str) -> bool:
    """
    Check if a file is in the ImmunoSEQ format based on its columns.

    :param file_path: Path to the file to check
    :return True if the file is an ImmunoSEQ file, False otherwise
    """
    df = pd.read_csv(file_path, sep='\t', nrows=0)  # Read only the header
    immunoseq_columns = [
        "nucleotide", 
        "aminoAcid", 
        "vGeneName", 
        "dGeneName", 
        "jGeneName", 
        "count (templates)", 
        "frequencyCount (%)"
    ]
    return any(col in df.columns for col in immunoseq_columns)


def parse_immunoseq(file_path: str) -> pd.DataFrame:
    """
    Parse an ImmunoSEQ file and rename its columns to the expected format.

    :param file_path: Path to the ImmunoSEQ file
    :return df: Dataframe with renamed columns
    """
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={
        'nucleotide': 'CDR3_nt',
        'aminoAcid': 'CDR3_aa',
        'vGeneName': 'V',
        'dGeneName': 'D',
        'jGeneName': 'J',
        'cGeneName': 'C',
        'count (templates/reads)': 'count',
        'frequencyCount (%)': 'proportion'
    })
    print(df.columns)
    return df


def parse(folder_path: str) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Load immune repertoire files from a folder and organize them into a 
    dictionary of dataframes.

    :param folder_path: Path to the folder containing the .tsv files
    :return: Tuple containing a dictionary of dataframes and a metadata dataframe
    """
    data: Dict[str, pd.DataFrame] = {}
    meta = []

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".tsv"):
            file_path = os.path.join(folder_path, file_name)
            
            # Determine if the file is in ImmunoSEQ format
            if is_immunoseq(file_path):
                df = parse_immunoseq(file_path)
            else:
                df = pd.read_csv(file_path, sep='\t')
            # Ensure all expected columns are present in the dataframe
            key_cols = [
                'CDR3_nt', 
                'CDR3_aa', 
                'V', 
                'D', 
                'J', 
                'C', 
                'count', 
                'proportion'
            ]

            for col in key_cols:
                if col not in df.columns:
                    df[col] = None  # Add missing columns with a default value of None
            
            for col in df.columns:
                if col not in key_cols:
                    df.drop(columns=[col], inplace=True)
            
            data[file_name] = df
            # Create metadata information for the file and add it to the meta list
            meta_info = {'Sample': file_name, 'Age': None, 'Status': None}
            meta.append(meta_info)
    
    # Convert the meta list to a pandas dataframe
    meta_df = pd.DataFrame(meta)
    return data, meta_df

# Usage example with dummy data
folder_path = '../pbmc-v2'
data, meta_df = parse(folder_path)


# Print out the data dictionary and meta dataframe
#print(data)
#print(list(data.keys()))
#print(data[list(data.keys())[1]])

