from preprocess import format_1, format_2
BIOKEY = format_1('raw_data/1882-BIOKEY_clonotypes_combined_cohort2.csv')
KCL_bulk, metadata = format_2('raw_data/pbmc-v2/')

from functions import plot_tcell_clonotype_totals, bar_plot_clonotype_proportions, do_all_bubbles




# ------ BIOKEY workflow ------------------------------------------

# plot_tcell_clonotype_totals(BIOKEY, 'SAMPLE_ID', 'figures/BIOKEY/clonotype_count/')
# bar_plot_clonotype_proportions(BIOKEY, 'SAMPLE_ID', 'figures/BIOKEY/clonotype_proportions/')
# pie_chart(BIOKEY, 'SAMPLE_ID', 'figures/BIOKEY/clonotype_proportions/')




# # ------ KCL_bulk workflow ----------------------------------------

# plot_tcell_clonotype_totals(KCL_bulk, 'patient', 'KCL-content/figures/clonotype_count/')
# bar_plot_clonotype_proportions(KCL_bulk, 'patient', 'KCL-content/figures/clonotype_proportions/')
# # # pie_chart(KCL_bulk, 'patient', 'KCL-content/figures/clonotype_proportions/')
do_all_bubbles()
