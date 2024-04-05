# regulonML
Sequence-based ML models of gene regulation

# instructions for including a new iModulon
1. See if iModulon is in data/TF_flags_expanded.csv, if not, add it
2. Check if genes are in data/saved_flags_expanded.csv, if not, add them
3. Need to fill out data/saved_flags_expanded.csv, steps listed:
    a. parameter_optimization/3_select_basal_condition.ipynb in order to pick basal conditions for the genes (should be one condition for all genes of an iModulon)
        i. Plug in [] to basal_conditions column to start
    b. 