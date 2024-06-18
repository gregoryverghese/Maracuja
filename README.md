<div style="text-align: left;">
<img src="maracuja.webp" alt="Maracuja Logo" width="100"/>
</div>

# Maracuja

Maracuja is a comprehensive package for parsing and analyzing TCR sequencing data. It supports data from various platforms, providing tools for detailed clonality analysis.

## Parsing

Maracuja supports parsing raw TCR sequencing files from Adaptive Biotechnologies and ImmunoSEQ. 

#### Supported Formats
- ImmunoSEQ

To parse the raw files, use the `rep_load` function:

```python
folder_path = 'path/to/your/folder_with_repertoire_files/'
data, meta_df = parser(folder_path)
```

## Analysis
Maracuja provides powerful tools to analyze clonality within TCR sequencing data.

### Clonality
You can visualize the top proportions of clones or rarest clones within your samples and stratify by a specific group

```python
top_prop = top_proportion(data)
plot_top_proportion_stacked(top_prop)

rare_prop = rare_proportion(data)
plot_top_proportion_grouped(rare_prop, meta_df, 'Group_Column_Name')
```

## License
Maracuja is licensed under the MIT License.
