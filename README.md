# Analysis

This repo/package contains the code which generates the content to store in KCL-Content which is then displayed in the TCR application here.

Package dependencies: matplotlib, pandas, numpy, circlify, sqlalchemy, palmotif, statistics, os, skunk, sklearn, re, pwseqdist, tcrdist. 

The general workflow to use the package is exemplified briefly in main.py. Comment and uncomment functions you require and run with `python main.py` 

### Initialise Database

On first clone, you must initialise the database using:

```python
initialise_db(<relative-path-to-bulkRNAseq-data>)
initialise_SC_db(<relative-path-to-scRNAseq-data>)
```

This saves each sample file from those directories as objects with attributes: id, pool, protocol, patient_id, HLA_A, and a function to retrieve data get_df() with columns: a_aminoacid (if applicable), b_aminoacid, count, proportion, and pMTnet (best prediction for all peptides from pool for both HLA-A types).

All sample objects can be accessed with `session.query(Sample).all()` and all objects are stored in a file generated called analysis.db.

### Bubble Overlay Plots

To create bubble overlay plots: use bubble_overlay() which requires the arguments: list of sample objects you wish to plot, list of x-axis order by patient and protocol,  list of y-axis order by pool, and the factor you wish to colour by including: none (no colour), EI (expansion index), TCT (total cell count), in_control (whether clonotype found in control) and pMTnet (by predicted affinity to peptides in pool). You can then either view plot with `plt.show()`, or save immediately with `plt.savefig(<path-to-save>, bbox_inches='tight')`

### Motif Analysis

Follow instructions on TCRdist website to install. If this works for you you can simply run in [main.py](http://main.py):

```python
arrays_folder = '../KCL-content/arrays'
logos_folder = '../KCL-content/logos'
plot_folder = '../KCL-content/figures/motif_analysis'
motifs1.aa_reads_CTR_array(samples, arrays_folder)
motifs2.arrays_to_logos_refs(arrays_folder, logos_folder) # To run in docker container if Palmotif is not able to install on your local
motifs3.plot_motifs('../KCL-content/refs_logos_each', x_order, y_order, plot_folder)
```

If you have problems building parasail, only the functions from motifs1 and motifs3 will work. TCRdist recommend using their docker image. 

```bash
docker run --name tcrdist -d -it -v $(pwd)/:/<dir-containing-both-Analysis-and-data-folder>/ quay.io/kmayerb/tcrdist3:0.1.9
docker exec -it tcrdist bash
```

### TCRdist

This part needs to be updated, but [TCRdist.py](http://TCRdist.py) only works if TCRdist install for you. TCR_MDS must be run outside of docker container.

### Notes

I am aiming eventually to change this workflow to not require the use of docker containers, or perhaps to use docker compose instead and have everything running in docker containers, so that everything can run cohesively.