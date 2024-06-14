import pandas as pd
from typing import Dict, List
from parser import parse


def clonal_proportion(data: Dict[str, pd.DataFrame], perc: int = 10) -> Dict[str, pd.Series]:
    """
    Calculate clonal proportion for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param perc: Percentage threshold to calculate clonal proportion.
    :return: Dictionary with filenames as keys and Series with clonal proportion data.
    """
    #print('greg')
    result = {}
    for sample, df in data.items():
        df = df.sort_values('count', ascending=False)
        total_count = df['count'].sum()
        cumulative_count = 0
        clone_count = 0
        #print(df)

        for count in df['count']:
            cumulative_count += count
            clone_count += 1
            if cumulative_count >= total_count * (perc / 100):
                break

        result[sample] = pd.Series({
            'Clones': clone_count,
            'Percentage': 100 * (cumulative_count / total_count),
            'Clonal.count.prop': clone_count / len(df)
        })
    return result

def clonal_space_homeostasis(data: Dict[str, pd.DataFrame], clone_types: Dict[str, float] = None) -> Dict[str, pd.Series]:
    """
    Calculate clonal space homeostasis for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param clone_types: Dictionary with clone type names and thresholds.
    :return: Dictionary with filenames as keys and Series with clonal space homeostasis data.
    """
    if clone_types is None:
        clone_types = {
            'Rare': 0.00001,
            'Small': 0.0001,
            'Medium': 0.001,
            'Large': 0.01,
            'Hyperexpanded': 1
        }
    clone_types = {'None': 0, **clone_types}
    
    result = {}
    for sample, df in data.items():
        proportions = df['count'] / df['count'].sum()
        homeostasis = {}
        
        for i in range(1, len(clone_types)):
            lower_bound = list(clone_types.values())[i - 1]
            upper_bound = list(clone_types.values())[i]
            homeostasis[f"{list(clone_types.keys())[i]} ({lower_bound} < X <= {upper_bound})"] = \
                proportions[(proportions > lower_bound) & (proportions <= upper_bound)].sum()
        
        result[sample] = pd.Series(homeostasis)
    return result


def top_proportion(data: Dict[str, pd.DataFrame], head: List[int] = [10, 100, 1000, 3000, 10000, 30000, 100000]) -> Dict[str, pd.Series]:
    """
    Calculate top clonal proportions for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param head: List of thresholds for top proportions.
    :return: Dictionary with filenames as keys and Series with top clonal proportions data.
    """
    result = {}
    for sample, df in data.items():
        df = df.sort_values('count', ascending=False)
        total_count = df['count'].sum()
        top_props = {}
        
        for h in head:
            top_props[f'Top {h}'] = df['count'].head(h).sum() / total_count
        
        result[sample] = pd.Series(top_props)
    return result


def rare_proportion(data: Dict[str, pd.DataFrame], bound: List[int] = [1, 3, 10, 30, 100]) -> Dict[str, pd.Series]:
    """
    Calculate rare clonal proportions for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param bound: List of thresholds for rare proportions.
    :return: Dictionary with filenames as keys and Series with rare clonal proportions data.
    """
    result = {}
    for sample, df in data.items():
        total_count = df['count'].sum()
        rare_props = {}
        
        for b in bound:
            rare_props[f'<= {b}'] = df['count'][df['count'] <= b].sum() / total_count
        rare_props['MAX'] = df['count'][df['count'] <= max(bound)].sum() / total_count
        
        result[sample] = pd.Series(rare_props)
    return result

# Usage example with dummy data

folder_path = '../pbmc-v2'
data, meta_df = parse(folder_path)

# Example usage of the functions
#clonal_prop = clonal_proportion(data, perc=10)
#clonal_homeo = clonal_space_homeostasis(data)
top_prop = top_proportion(data)
#rare_prop = rare_proportion(data)

# Print results
#print(clonal_prop)
#print(clonal_homeo)
print(top_prop)
#print(rare_prop)

