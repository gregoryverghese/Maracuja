# Analysis

This repo/package contains the code which generates the content to store in KCL-Content which is then displayed in the TCR application.

Package dependencies: matplotlib, pandas, numpy, circlify, sqlalchemy, palmotif, statistics, os, skunk, sklearn, re, pwseqdist, tcrdist. 

### Process Raw data

Process feature count data into standardised form with can be added as object to database. Run process_raw_data.py followed by the format of the data, 1 for adaptive data, 2 for data from immunoseq, and the path to the data files. You must manually create the metadata file describing the protocol, patient info, and type of data.

### Initialise Database

On first clone, you must initialise the database using:

```python
initialise_db(<relative-path-to-bulkDNAseq-data>)
initialise_SC_db(<relative-path-to-scRNAseq-data>)
```

This saves each sample file from those directories as objects with attributes: id, pool, protocol, patient_id, HLA_A, and a function to retrieve data get_df() with columns: a_aminoacid (if applicable), b_aminoacid, count, proportion, and pMTnet (best prediction for all peptides from pool for both HLA-A types).

All sample objects can be accessed with `session.query(Sample).all()` and all objects are stored in a file generated called analysis.db.

### Bubble Plots

To create bubble plots for each sample, run bubbles.py followed by path and outdir to generate bubble plots for every sample in a given directory.

To create bubble overlay plots: use bubble_overlay() which requires the arguments: list of sample objects you wish to plot, list of x-axis order by patient and protocol,  list of y-axis order by pool, and the factor you wish to colour by including: none (no colour), EI (expansion index), TCT (total cell count), in_control (whether clonotype found in control) and pMTnet (by predicted affinity to peptides in pool). You can then either view plot with `plt.show()`, or save immediately with `plt.savefig(<path-to-save>, bbox_inches='tight')`

### Summary stats

Run summary_stats.py followed by the path to data and the path to KCL-Content/summary_stats to generate the top 5 clonotypes plus their comparative proportions in the same sample controls.

### Shannon's index

Run shannons_index.py followed by the path to data to calculate Shannon's index and simpson's index and to statistically verify this between each treatment sample and its control.

