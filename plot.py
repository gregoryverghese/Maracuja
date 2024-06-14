import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Dict

from clonality import *

def plot_top_proportion(result: Dict[str, pd.Series]):
    """
    Plot the top clonal proportions for each sample.

    :param result: Dictionary with filenames as keys and Series with top clonal proportions data.
    """
    # Convert the result dictionary to a DataFrame for plotting
    df = pd.DataFrame(result).T.reset_index().rename(columns={'index': 'Sample'})
    
    # Melt the DataFrame to long format for Seaborn plotting
    df_melted = df.melt(id_vars=['Sample'], var_name='Top Clones', value_name='Proportion')
    
    # Plotting
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Top Clones', y='Proportion', hue='Sample', data=df_melted)
    plt.title('Top Clonal Proportions')
    plt.xlabel('Top Clones')
    plt.ylabel('Proportion')
    plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Usage example
folder_path = '../pbmc-v2'
data, meta_df = parse(folder_path)

top_prop = top_proportion(data)
plot_top_proportion(top_prop)

