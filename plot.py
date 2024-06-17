import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Dict


def plot_top_proportion_grouped(result: Dict[str, pd.Series], metadata: pd.DataFrame, group_by: str):
    """
    Plot the top clonal proportions for each sample as a grouped bar chart.

    :param result: Dictionary with filenames as keys and Series with top clonal proportions data.
    :param metadata: DataFrame with metadata information for the samples.
    :param group_by: Column in the metadata to use for grouping.
    """
    # Convert the result dictionary to a DataFrame for plotting
    df = pd.DataFrame(result).T.reset_index().rename(columns={'index': 'Sample'})
    df_melted = df.melt(id_vars=['Sample'], var_name='Top Clones', value_name='Proportion')
    # Merge with metadata to get grouping information
    df_melted = df_melted.merge(metadata, left_on='Sample', right_on='Sample')
  
    
    # Melt the DataFrame to long format for better plotting
    
    #print(df_melted)
    
    # Plotting
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Top Clones', y='Proportion', hue=group_by, data=df_melted, ci=None)
    
    # Set plot title and labels
    plt.title('Top Clonal Proportions (Grouped)')
    plt.xlabel('Top Clones')
    plt.ylabel('Proportion (%)')
    
    # Place legend outside of the plot
    plt.legend(title=group_by, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Ensure tight layout
    plt.tight_layout()
    
    # Show plot
    plt.show()
    
    return df_melted


def plot_top_proportion_stacked(result: Dict[str, pd.Series]):
    """
    Plot the top clonal proportions for each sample as a percent-stacked bar plot.

    :param result: Dictionary with filenames as keys and Series with top clonal proportions data.
    """
    # Convert the result dictionary to a DataFrame for plotting
    df = pd.DataFrame(result).T.reset_index().rename(columns={'index': 'Sample'})

    #df_normalized = df * 100

    # Plotting
    ax = df.set_index('Sample').plot(kind='bar', stacked=True, figsize=(12, 8), colormap=cm.get_cmap('viridis'))

    # Set plot title and labels
    plt.title('Top Clonal Proportions (Percent-Stacked)')
    plt.xlabel('Sample')
    plt.ylabel('Proportion (%)')

    # Place legend outside of the plot
    plt.legend(title='Top Clones', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Ensure tight layout
    plt.tight_layout()

    # Show plot
    plt.show()
