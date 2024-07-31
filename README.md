# regulonML
Sequence-based ML models of gene regulation

# How to get started
1. Install necessary packages - XXX figure out a way to list them all
*Note - All files are already included for utilizing P1K (E. coli iModulons). If wanting to use these, skip to step XXX.
2. Download external data files, move to data/external/. These include:
    a. data/imodulon_info/ - Needs to contain A.csv, M.csv, iM_table.csv, log_tpm.csv. Optionally can include sample_table.csv, gene_info.csv for use in some troubleshooting. Necessary files can most easily be found on https://imodulondb.org.
    b. data/regulonDB_files - Needs to contain network_tf_operon.txt, NetWorkTFGene.txt, OperonSet.txt. These can be downloaded from https://regulondb.ccg.unam.mx/datasets. XXX create functionality to create TRN_df.pkl from these inputs (steal code from TF KALE I believe).
    c. data/validation_data_sets/ - Various data sets used to validate models. Only heinemann_protein_conc.xlsx is necessary, which provides limitations on model parameters and is available in the supplement of 10.1038/nbt.3418. Optionally can include variuos other validation datasets, most importantly metabolite data which by default is in the "stationary_phase" folder sourced from 10.1038/nmeth.3584.
    d. (optional) data/external/KEGG/ - Need pathway_to_gene.txt and pathways.txt, used in some troubleshooting
3. Set up the option files in data/options/
    a. data/options contains two files:
        i. gene_flags.csv contains options for how the model handles each of the individual genes.
        ii. TF_flags.csv contains options for how the model handles each individual regulator group.
    b. First, run 1_model_setup/1_setup_data_options.
        i. This will generate semi-complete flags files in data/options titled template_gene_flags.csv and template_TF_flags.csv
    c. Set the effector and effector type for each regulator in template_TF_flags.csv that you intend to include:
        i. (note) - Both a small test case and more complete case of regulators for P1K (E. coli iModulons) is included already, so no need to complete this manual step if you plan on using P1K.
        ii. Set the effector, make sure the naming matches that of any included metabolites data. By default, these names are in validation_data_sets/stationary_phase/metaboltie_data.xlsx.
        iii. The following options are currently available, for each regutor modify one column that best fits it to 1:
            1. cAct_no_effector_form - TF is a promoter which activates without needing to bind an effector
            2. cAct_multi_effector_binding - TF is a promoter which activates with the binding of one or more copies of the same effector
            3. cAct_multi_co_effector_binding - TF is a promoter which activates with the binding of two different effectors
            4. cInh_no_effector_form - TF is an inhibitor which activates without needing to bind an effector
            5. cInh_multi_effector_binding - TF is an inhibitor which activates with the binding of one or more copies of the same effector
            6. cInh_multi_co_effector_binding -  TF is an inhibitor which activates with the binding of two different effectors
    d. Rename the template files to gene_flags.csv and TF_flags.csv and place them in the data/options folder.
4. Set basal conditions for the individual genes
    a. Run 1_model_setup/2_select_basal_conditions.ipynb
        i. This will set the same basal condition for all genes of a regulatory case. The basal condition should correspond to zero regulator activity. Gene expression is standardized for each gene across all samples. The basal conditions is selected based on if the case's regulator is a promoter, inhibitor, or dual-regulator:
            1. Regulator is a promoter - The sample with the lowest expression is picked
            2. Regulator is an inhibitor - The sample with the highest expression is picked
            3. Regulator is a dual-regulator - The sample with the most average expression is picked
        ii. You may want to manually set some basal conditions, if you have samples that are better approximations of zero activity, such as knockouts. Basal conditions need to be the same across all genes of the same regulatory case.
5. Select grid values for each gene.
    a. Run 1_model_setup/3_pick_grid_for_gene.ipynb
        i. Each gene needs a set of associated constants (KdRNAP, KeqOpening, KdRNAPActivator) to solve the model. This notebook generates these constants for each.
        ii. 9 different sets of solutions are generated ranging from low to high KdRNAP values. Typically grid #7 works best and is the default.
        iii. KdRNAPActivator is needed for any case with an activator. The model can pick an ideal KdRNAPActivator value if one exists, but this process often has no solution and occasionally fails. This requires some manual work, so for each gene check the following:
            1. The output for each gene will include the gene name followed by its activating and inhibiting iModulon.
            2. Check if the range for mRNA ratio makes sense.
                a. mRNA ratio should be almost entirely below 1 if inhibitor only, almost entirely above 1 if promoter only, mean should be near 1 if dual-regulator.
            3. If inhibitor - Check for negative correlation between cInhibitor and actual mRNA ratio, normal range for cInhibitor values (0-1000 is ideal, 0-10000 is okay, 0-100000 may cause overfitting on outliers in later steps).
            4. If activator - Need to verify the KdRNAPActivator value:
                a. 1st Objecitive Function and 2nd Objective function are used to pick a KdRNAPActivator value which balances between a KdRNAPActivator value which causes the activator increase recruitment of RNAP to the promoter site while not making it a completely binary switch.
                b. Ideally, the vertical dotted line will be at the minimum point in both plots. If it is not, there are a few options:
                    i. Look at the same plots for all 9 grid values, if it works for one of these, select this grid value.
                    ii. If the dotted line is not picking the minimum point but there appears to be one for one of the grid values - modify the initial_guess_ratio column in data/options/gene_flags.csv to somewhere between 0 and 1 and rerun the code to regenerate the results. This sets an initial different guess for KdRNAPActivator which often allows the optimization code to work properly in this case.
            5. If activator and valid KdRNAPActivator value - Check for positive correlation between cActivator and actual mRNA ratio, normal range for cActivator values (0-1000 is ideal, 0-10000 is okay, 0-100000 may cause overfitting on outliers in later steps).
            6. If previous steps worked, change the "checked" and "include" columns in the gene_flags.csv to TRUE. Set the "grid_use" column to the preferred grid value.
        iv. Generally, for the above steps, I found it easier to iterate over everything once and then reiterate through the genes which failed the checks. If a gene is unable to work, you may simply change the "checked" column to TRUE and the "include" column to FALSE. This will mean the gene will not be included in the model.
6. Generate necesarry files for model generation.
    a. Run 1_model_setup/4_generate_model_files.ipynb, which will iterate over each gene and generate the necessary files.
7. Optimize the GAMS model.
    a. Run 2_GAMS/1_run_model.ipynb, which will run GAMS to optimize the model for the input data.
8. Visualize results.
    a. Run various notebooks in 3_visualize_results to look at the GAMS results.
9. Optimize parameters if necessary.
    a. Run 4_optimization/1_optimize_GAMS_parameters.ipynb to select better constraints on the parameters for GAMS to create more accurate models.
    b. This process is very slow and is best run either during computer downtime or on an independent server. It involves running GAMS on a wide variety of parameter constraints, taking the best result, creating a new set of parameter results closer to the best result, and repeating this process a set number of times.




# Terminology
1. Regulatory case - A set of genes which all have the same promoter and inhibitor. For example: All genes inhibited by argR and not promoted by any other genes are a regulatory case.
2. cInhibitor - XXX
3. cActivator - XXX



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