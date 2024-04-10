# regulonML
Sequence-based ML models of gene regulation

# instructions for including a new iModulon
1. See if iModulon is in data/TF_flags_expanded.csv, if not, add it, fill out, copy over to data/TF_saved_flags.csv
2. Check if genes are in data/saved_flags_expanded.csv, if not, add them
3. Need to fill out data/saved_flags_expanded.csv, steps listed:
    a. parameter_optimization/3_select_basal_condition.ipynb in order to pick basal conditions for the genes (should be one condition for all genes of an iModulon)
        i. If iM is activator, pick lowest condition, if inhibitor, pick highest
        ii. Take saved off data/exported_flags_df.csv, add it to bottom of flags_df.csv
    b. Copy down the rest of the missing columns
    c. Reupload saved_flags.csv
4. Need to pick grid conditions in parameter_optimization/2_pick_grid_for_gene.ipynb
    a. modify iM in first cell
    b. set gene_int to 0 in second cell, increase by one each iteration
    c. run and look at which grid_ct looks best, set as such in the saved_flags.csv
        i. it's almost always 7, sometimes 8
    d. repeat for all genes
    e. Reupload saved_flags.csv