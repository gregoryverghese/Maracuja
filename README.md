
# Maracuja ![Maracuja Logo](maracuja.png)

Maracuja is a comprehensive package for parsing and analyzing TCR sequencing data. It supports data from various platforms, providing tools for detailed clonality analysis.

## Table of Contents
- [Parsing](#parsing)
- [Analysis](#analysis)
  - [Clonality](#clonality)
    - [Stacked Bar Chart](#stacked-bar-chart)
    - [Grouped Bar Chart](#grouped-bar-chart)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Parsing

Maracuja supports parsing raw TCR sequencing files from Adaptive Biotechnologies and ImmunoSEQ. 

### Supported Formats
- **Adaptive Biotechnologies**
- **ImmunoSEQ**

### Parsing Example
To parse the raw files, use the `rep_load` function:

```python
folder_path = 'path/to/your/folder_with_repertoire_files/'
data, meta_df = rep_load(folder_path)

## Analysis
Maracuja provides powerful tools to analyze clonality within TCR sequencing data.

### Clonality
You can create a stacked bar chart to visualize the top proportions of clones within your samples.

```python
top_prop = top_proportion(data)
plot_top_proportion_stacked(top_prop)

Grouped Bar Chart
Stratify your samples based on a particular group to create a grouped bar chart using rare proportions.

```python
rare_prop = rare_proportion(data)
plot_top_proportion_grouped(rare_prop, meta_df, 'Group_Column_Name')

## License
Maracuja is licensed under the MIT License.
