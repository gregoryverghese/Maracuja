import os

import pandas as pd
import random
from typing import Dict, Tuple,  List
import matplotlib.pyplot as plt
import seaborn as sns


def top_proportion(data: Dict[str, pd.DataFrame], head: List[int] = [10, 100, 1000, 3000, 10000, 30000, 100000]) -> Dict[str, pd.Series]:
    """
    Calculate absolute top clonal proportions for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param head: List of thresholds for top proportions.
    :return: Dictionary with filenames as keys and Series with absolute top clonal proportions data.
    """
    result = {}
    for sample, df in data.items():
        df = df.sort_values('count', ascending=False)
        total_count = df['count'].sum()
        top_props = {}

        previous_threshold = 0
        for h in head:
            if previous_threshold < len(df):
                current_threshold = min(h, len(df))
                top_props[f'Top {h}'] = df['count'].iloc[previous_threshold:current_threshold].sum() / total_count
                previous_threshold = current_threshold
            else:
                top_props[f'Top {h}'] = 0

        result[sample] = pd.Series(top_props)
    return result


def rare_proportion(data: Dict[str, pd.DataFrame], bound: List[int] = [1, 3, 10, 30, 100]) -> Dict[str, pd.Series]:
    """
    Calculate absolute rare clonal proportions for each sample in the data.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param bound: List of thresholds for rare proportions.
    :return: Dictionary with filenames as keys and Series with absolute rare clonal proportions data.
    """
    result = {}
    for sample, df in data.items():
        total_count = df['count'].sum()
        rare_props = {}

        previous_threshold = 0
        for b in bound:
            current_threshold = b
            if previous_threshold < len(df):
                rare_props[f'<= {b}'] = df['count'][(df['count'] > previous_threshold) & (df['count'] <= current_threshold)].sum() / total_count
            else:
                rare_props[f'<= {b}'] = 0
            previous_threshold = current_threshold

        # Adding a category for the maximum bound
        if max(bound) < len(df):
            rare_props[f'>{max(bound)}'] = df['count'][df['count'] > max(bound)].sum() / total_count
        else:
            rare_props[f'>{max(bound)}'] = 0

        result[sample] = pd.Series(rare_props)
    return result
