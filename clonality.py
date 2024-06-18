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


def bubble_overlay(data: Dict[str, pd.DataFrame], x_order: List[str], y_order: List[str], factor: str, metadata: pd.DataFrame):
    """
    Create a grid of bubble plots for multiple samples.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param x_order: List defining the order of x-axis categories.
    :param y_order: List defining the order of y-axis categories.
    :param factor: Factor determining the color of the bubbles.
    :param metadata: DataFrame with metadata information for the samples.
    """

    for sample_name, df in data.items():
        try:
            sample_metadata = metadata[metadata['Sample'] == sample_name].iloc[0]
        except:
            continue
        patient_protocol = f"{sample_metadata['patient_id']}"
        
        x_row = x_order.index(sample_metadata['pool_id'])
        y_col = y_order.index(patient_protocol)

        ax = plt.subplot2grid((len(y_order),len(x_order)), (y_col,x_row))
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
        ax.axis('off')

        df = data[sample_name].head(100)
        expansion_index = df['proportion'].max()
        colour = [(expansion_index, 1 - expansion_index, 0)] * len(df)
       
        # Prepare data for circlify
        data_prepared = [{'id': str(i), 'datum': df['count'][i]} for i in df.index]
        circles = circ.circlify(data_prepared)
        print('circle done')

        # Add circles to the plot
        for circle, col in zip(circles, df['col']):
            c = Circle((circle.x, circle.y), circle.r, color=col, linewidth=0)
            ax.add_patch(c)

    ax = plt.subplot2grid((len(y_order),len(x_order)), (0,0))
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')

    for i in range(len(y_order)):
        plt.text(-6, -2.05*i, y_order[i], fontsize = 10)
    for i in range(len(x_order)):
        plt.text(2.48*i-0.5, -35, x_order[i], fontsize = 10, rotation=90)

    plt.text(9, -37, 'Treatment', fontsize = 10)
    plt.text(-7, -8, 'Patient', fontsize = 10, rotation=90)
    ax.axis('off')



def bubble_overlay(data: Dict[str, pd.DataFrame], x_order: List[str], y_order: List[str], factor: str, metadata: pd.DataFrame):
    """
    Create a grid of bubble plots for multiple samples.

    :param data: Dictionary with filenames as keys and dataframes as values.
    :param x_order: List defining the order of x-axis categories.
    :param y_order: List defining the order of y-axis categories.
    :param factor: Factor determining the color of the bubbles.
    :param metadata: DataFrame with metadata information for the samples.
    """

    for sample_name, df in data.items():
        print(sample_name)
        print(metadata[metadata['Sample'] == sample_name])
        try:
            sample_metadata = metadata[metadata['Sample'] == sample_name].iloc[0]
        except:
            continue
        patient_protocol = f"{sample_metadata['patient_id']}"
        
        x_row = x_order.index(sample_metadata['pool_id'])
        y_col = y_order.index(patient_protocol)

        ax = plt.subplot2grid((len(y_order),len(x_order)), (y_col,x_row))
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')
        ax.axis('off')

        df = data[sample_name].head(100)
        expansion_index = df['proportion'].max()
        colour = [(expansion_index, 1 - expansion_index, 0)] * len(df)
       
        # Prepare data for circlify
        data_prepared = [{'id': str(i), 'datum': df['count'][i]} for i in df.index]
        circles = circ.circlify(data_prepared)
        print('circle done')

        # Add circles to the plot
        for circle, col in zip(circles, df['col']):
            c = Circle((circle.x, circle.y), circle.r, color=col, linewidth=0)
            ax.add_patch(c)

    ax = plt.subplot2grid((len(y_order),len(x_order)), (0,0))
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal')

    for i in range(len(y_order)):
        plt.text(-6, -2.05*i, y_order[i], fontsize = 10)
    for i in range(len(x_order)):
        plt.text(2.48*i-0.5, -35, x_order[i], fontsize = 10, rotation=90)

    plt.text(9, -37, 'Treatment', fontsize = 10)
    plt.text(-7, -8, 'Patient', fontsize = 10, rotation=90)
    ax.axis('off')


