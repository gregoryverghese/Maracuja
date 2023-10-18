import matplotlib.pyplot as plt
from database.initialise_db import Sample, SC_Sample, session, initialise_db, initialise_SC_db
from bubbles.bubbles import plot
import motifs.motifs1 as motifs1
### import motifs.motifs2 as motifs2 # modules used here doen't work on my local, so I use the docker image. 
import motifs.motifs3 as motifs3


### Initialise Data
### Only needs to be done once, or again when data is added/change
initialise_db('../KCL-Raw-data')
initialise_SC_db('../Demo-SC-data')



### BUBBLE PLOT ALL
samples = [sample for sample in session.query(Sample).all() if sample.protocol != 'T' ]
x_order = ['CEF', 'DMSO', 'Neg', 'Pool1', 'Pool2', 'Pool5', 'Pool7', 'Pool8', 'Pool9', 'Pool10']
y_order = ['KCL710-A', 'KCL717-A', 'KCL725-A', 'KCL737-A', 'KCL750-A', 'KCL763-A', 'KCL710-B', 'KCL763-B']


### BUBBLE PLOT SC
samples_SC = [sample for sample in session.query(SC_Sample).all() if sample.protocol != 'T' ]
x_order_SC = ['Pool1', 'Pool2']
y_order_SC = ['KCL580-B', 'KCL500-B']


### PLOT
# bubble_overlay(samples_SC, x_order_SC, y_order_SC, 'none')
# bubble_overlay(samples, x_order, y_order, 'in_control')
# plt.show()

### MOTIF ANALYSIS
arrays_folder = '../KCL-content/arrays'
logos_folder = '../KCL-content/logos'
plot_folder = '../KCL-content/figures/motif_analysis'
# motifs1.aa_reads_CTR_array(samples, arrays_folder)
# ### motifs2.arrays_to_logos_refs(arrays_folder, logos_folder) # To run in docker container if Palmotif is not able to install on your local
# motifs3.plot_motifs('../KCL-content/refs_logos_each', x_order, y_order, plot_folder)