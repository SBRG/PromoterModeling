"""
Contains functions pertaining to the running and reading of GAMS

Functions:
run_multi_GAMS - Runs GAMS on all genes at once of a specific regulatory case
read_multi_GAMs - Reads GAMS results on all genes at once of a specific regulatory case
"""

# imports
import subprocess
import os
import pandas as pd
import numpy as np
import pickle
import shutil
import platform
import itertools


# something that likely needs to be changed if moving away from E. coli
iM_to_b_regulator = {
    'Purine' : ['b1658'], 
    'Sugar Diacid': ['b0162'],
    'Translation': ['b0145'],
    'OxyR': ['b3961'],
    'FlhDC-2': ['b1892'],
    'Glutamine': ['b3868'],
    'RpoH': ['b3461'],
    'Methionine': ['b3938'],
    'SoxS': ['b4062'],
    'YjfJ': ['b4182'],
    'Lrp': ['b0889'],
    'CpxR': ['b3912'],
    'EvgA': ['b2369'],
    'Biotin': ['b3973'],
    'BasR': ['b4113'],
    'Sulfoquinovose': ['b3884'],
    'RpoS': ['b2741'],
    'Fur-1': ['b0683'],
    'Fur-2': ['b0683'],
    'Magnesium': ['b1130'],
    'Fnr-1': ['b1334'],
    'Zinc-2': ['b4004'],
    'NRZ': ['b3405'],
    'RpoE': ['b2573'],
    'Phosphate-1': ['b0399'],
    'ArcA': ['b4401'],
    'Potassium': ['b0694'],
    'DNA Damage': ['b4043'],
    'YieP': ['b3755'],
    'Nitrogen': ['b1988'],
    'Crp-1': ['b3357'],
    'YgeV': ['b2869'],
    'FDH-O': ['b4764'],
    'YgbI': ['b2735'],
    'Lysine/T2SS': ['b2916'],
    'Fnr-3': ['b1334'],
    'FliA': ['b1922'],
    'Pyruvate-1': ['b2125'],
    'Crp-2': ['b3357'],
    'Cysteine-1': ['b1275'],
    'Pyruvate-2': ['b2381'],
    'Phosphate-2': ['b0399'],
    'Capsule': ['b1951'],
    'FlhDC-1': ['b1892'],
    'Arginine': ['b3237'],
    'Nickel/Cobalt': ['b2105'],
    'YmfT': ['b1146'],
    'Fnr-2': ['b1334'],
    'Cysteine-2': ['b1275'],
    'YcjW': ['b1320'],
    'NrdR': ['b0413'],
    'Cra': ['b0080'],
    'Maltose': ['b3418'],
    'Fatty Acid': ['b1187'],
    'DhaR': ['b1201'],
    'gcvB': ['b4443']
}


def run_multi_GAMS(flags_df, TF_flags_df, stable_flags, cell_constants, GAMS_run_dir, GAMS_exec = 'gams', parameter_flags = None):
    """
    Inputs various flags and constants, runs GAMS, saves to directed folder.
    
    Inputs:
        flags_df (dict) : settings flags specific to the genes
        TF_flags_df (dict) : settings flags specific to the regulators
        stable_flags (dict) : dictionary of settings flags
        cell_constants (dict) : set of biological constants
        GAMS_run_dir (string) : output location for GAMS run
        GAMS_exec (string) : path to GAMS executable, if not set defaults to assume GAMS can be called from command line by "gams"
        parameter_flags (dict) : overwrites default parameters if set
    
    Returns:
        None
    """
    
    ############################################################
    # save input parameters
    ############################################################
    if parameter_flags is None:
        parameter_flags = stable_flags
        baby_dict = {
            'act_TF_conc_lo' : parameter_flags['act_TF_conc_lo'],
            'act_TF_conc_up' : parameter_flags['act_TF_conc_up'],
            'act_Kd_lo' : parameter_flags['act_Kd_lo'],
            'act_Kd_up' : parameter_flags['act_Kd_up'],
            'inh_TF_conc_lo' : parameter_flags['inh_TF_conc_lo'],
            'inh_TF_conc_up' : parameter_flags['inh_TF_conc_up'],
            'inh_Kd_lo' : parameter_flags['inh_Kd_lo'],
            'inh_Kd_up' : parameter_flags['inh_Kd_up'],

            'weight_act_obj1' : parameter_flags['weight_act_obj1'],
            'weight_inh_obj1' : parameter_flags['weight_inh_obj1'],
            'weight_act_obj2' : parameter_flags['weight_act_obj2'],
            'weight_inh_obj2' : parameter_flags['weight_inh_obj2'],
            'weight_mRNA_match' : parameter_flags['weight_mRNA_match'],
            'weight_act_corr' : parameter_flags['weight_act_corr'],
            'weight_inh_corr' : parameter_flags['weight_inh_corr'],
            'act_metab_Total_lo' : parameter_flags['act_metab_Total_lo'],
            'act_metab_Total_up' : parameter_flags['act_metab_Total_up'],
            'inh_metab_Total_lo' : parameter_flags['inh_metab_Total_lo'],
            'inh_metab_Total_up' : parameter_flags['inh_metab_Total_up'],
        }
    else:
        para_sweep = ['act_TF_conc_lo', 'act_TF_conc_up', 'act_Kd_lo', 'act_Kd_up', 'inh_TF_conc_lo', 'inh_TF_conc_up', 'inh_Kd_lo', 'inh_Kd_up', 'weight_act_obj1', 'weight_inh_obj1', 'weight_act_obj2', 'weight_inh_obj2', 'weight_mRNA_match', 'weight_act_corr', 'weight_inh_corr', 'inh_metab_Total_lo', 'inh_metab_Total_up', 'act_metab_Total_lo', 'act_metab_Total_up']
        baby_dict = {}
        for para in para_sweep:
            if para in parameter_flags:
                baby_dict.update({para : parameter_flags[para]})
    df = pd.DataFrame(list(baby_dict.items()), columns=['Parameter', 'Value']).set_index('Parameter')
    df.to_csv(GAMS_run_dir+'/input_files/parameters.csv')
    
    
    ############################################################
    # TF specific parameters
    ############################################################
    TF_flags_df.to_csv(GAMS_run_dir+'/input_files/TF_constants.csv')
    
    
    ############################################################
    # gather files needed to run GAMS
    ############################################################
    files_use = []
    for sample in stable_flags['limit_samples']:
        use_zerod_A = flags_df.loc[sample]['use_zerod_A_matrix']
        if stable_flags['use_greedy']:
            f_name = sample+'_greedy.pkl'
        else:
            f_name = sample+'.pkl'
        if os.path.exists('../data/processed/cAct_cInh_vals/'+f_name):
            files_use.append(f_name)
    
    for f in files_use:
        shared_indices = []
        first = True
        for file in files_use:
            gene_name = file.split('.')[0].split('_')[0]
            temp_df = pd.read_pickle('../data/processed/cAct_cInh_vals/'+file)
            if first:
                shared_indices = set(temp_df.index)
                first = False
            else:
                shared_indices = set(shared_indices.intersection(set(temp_df.index)))
        shared_indices = list(shared_indices)
        shared_indices.sort()
        act_df = pd.DataFrame(index = shared_indices)
        inh_df = pd.DataFrame(index = shared_indices)
        for file in files_use:
            gene_name = file.split('.')[0].split('_')[0]
            act_df[gene_name] = pd.read_pickle('../data/processed/cAct_cInh_vals/'+file)['cAct'].loc[act_df.index].values
            inh_df[gene_name] = pd.read_pickle('../data/processed/cAct_cInh_vals/'+file)['cInh'].loc[inh_df.index].values
    
    # pick samples to use
    samples_use = act_df.index.to_list()
    if 'small_dataset' in stable_flags and stable_flags['small_dataset'] == True:
        stationary_samples = ['starve_series__t00_growth1',
            'starve_series__t01_starve',
            'starve_series__t22_growth2',
            'starve_series__t23_growth2',
            'starve_series__t27_growth2',
            'starve_series__t28_growth2',
            'starve_series__t06_starve',
            'starve_series__t08_starve',
            'starve_series__t09_starve',
            'starve_series__t10_starve',
            'starve_series__t11_starve',
            'starve_series__t12_starve',
            'starve_series__t16_starve',
            'starve_series__t17_starve',
            'starve_series__t18_starve',]
    
        samples_use = samples_use[::-1][0:35]+stationary_samples
        samples_use = list(set(samples_use).intersection(act_df.index).intersection(inh_df.index))
    act_df = act_df.loc[samples_use]
    inh_df = inh_df.loc[samples_use]
    
    # save off max values for scaling objective functions
    vals = []
    for col in act_df.columns:
        vals.append(max(act_df[col]))
    max_df = pd.DataFrame(vals, index = act_df.columns)
    max_df.to_csv(GAMS_run_dir+'/input_files/max_cActs.csv')
    vals = []
    for col in inh_df.columns:
        vals.append(max(inh_df[col]))
    max_df = pd.DataFrame(vals, index = inh_df.columns)
    max_df.to_csv(GAMS_run_dir+'/input_files/max_cInhs.csv')
    
    # convert this into iModulon format that looks like excel so GAMS can read it
    # also need to create mapping between existing gene- iM in act and gene-iM in inh pairs
    act_gene_iM_pairs = []
    inh_gene_iM_pairs = []
    act_df = act_df.T
    inh_df = inh_df.T
    iMs = []
    keep = []
    for gene in act_df.index:
        iM = flags_df.loc[gene]['act_iM']
        if type(iM) == float:
            keep.append(False)
        else:
            act_gene_iM_pairs.append((gene, iM))
            keep.append(True)
            iMs.append(iM)
    act_df = act_df.loc[keep]
    cols = inh_df.columns.to_list()
    act_df['iM'] = iMs
    cols.insert(0, 'iM')
    act_df['iM'] = iMs
    act_df = act_df[cols]

    iMs = []
    keep = []
    for gene in inh_df.index:
        iM = flags_df.loc[gene]['inh_iM']
        if type(iM) == float:
            keep.append(False)
        else:
            inh_gene_iM_pairs.append((gene, iM))
            keep.append(True)
            iMs.append(iM)
    inh_df = inh_df.loc[keep]
    cols = inh_df.columns.to_list()
    cols.insert(0, 'iM')
    inh_df['iM'] = iMs
    inh_df = inh_df[cols]
    
    
    
    # include missing values as zeros and finish making mapping df
    all_iMs = list(set(act_df.iM).union(inh_df.iM))
    genes = list(set(act_df.index).union(set(inh_df.index)))
    needs_to_exist = [gene+';'+iM for gene, iM in itertools.product(genes, all_iMs)]
    
    # add missing to cAct
    exists = [gene+';'+iM for gene, iM in zip(act_df.index, act_df.iM)]
    to_create = list(set(needs_to_exist) - set(exists))
    new_row_indices = [val.split(';')[0] for val in to_create]
    new_iMs = [val.split(';')[1] for val in to_create]
    new_rows = pd.DataFrame(index = new_row_indices, columns = act_df.columns)
    new_rows['iM'] = new_iMs
    act_df = pd.concat([act_df, new_rows])
    act_df = act_df.fillna(0)

    # add missing to cAct
    exists = [gene+';'+iM for gene, iM in zip(inh_df.index, inh_df.iM)]
    to_create = list(set(needs_to_exist) - set(exists))
    new_row_indices = [val.split(';')[0] for val in to_create]
    new_iMs = [val.split(';')[1] for val in to_create]
    new_rows = pd.DataFrame(index = new_row_indices, columns = inh_df.columns)
    new_rows['iM'] = new_iMs
    inh_df = pd.concat([inh_df, new_rows])
    inh_df = inh_df.fillna(0)
    
    act_df.to_csv(GAMS_run_dir+'/input_files/composite_cAct_vals.csv')
    inh_df.to_csv(GAMS_run_dir+'/input_files/composite_cInh_vals.csv')
    
    # make mapping
    act_mapping_df = pd.DataFrame(index = list(set(act_df.index)), columns = all_iMs).fillna(0)
    inh_mapping_df = pd.DataFrame(index = list(set(inh_df.index)), columns = all_iMs).fillna(0)
    for gene, iM in act_gene_iM_pairs:
        act_mapping_df.at[gene, iM] = 1
    for gene, iM in inh_gene_iM_pairs:
        inh_mapping_df.at[gene, iM] = 1
    act_mapping_df.to_csv(GAMS_run_dir+'/input_files/cAct_mapping.csv')
    inh_mapping_df.to_csv(GAMS_run_dir+'/input_files/cInh_mapping.csv')
    
    # if the legnth of iMs is 1, we need to add a comma to get it to run here
    if len(all_iMs) == 1:
        with open(GAMS_run_dir+'/input_files/cAct_mapping.csv', 'r') as f:
            lines = f.readlines()
        lines[0] = lines[0].rstrip('\n')+','+'\n'
        with open(GAMS_run_dir+'/input_files/cAct_mapping.csv', 'w') as f:
            f.writelines(lines)
        with open(GAMS_run_dir+'/input_files/cInh_mapping.csv', 'r') as f:
            lines = f.readlines()
        lines[0] = lines[0].rstrip('\n')+','+'\n'
        with open(GAMS_run_dir+'/input_files/cInh_mapping.csv', 'w') as f:
            f.writelines(lines)
    
    # need to save off dummy dimensional dataframe
    all_samples = list(set(act_df.columns).union(inh_df.columns))
    all_samples.remove('iM')
    all_samples.insert(0, 'iM')
    new_row_indices = [val.split(';')[0] for val in needs_to_exist]
    new_iMs = [val.split(';')[1] for val in needs_to_exist]
    dummy_df = pd.DataFrame(index = new_row_indices, columns = all_samples)
    dummy_df['iM'] = new_iMs
    dummy_df = dummy_df.fillna(1)
    dummy_df.to_csv(GAMS_run_dir+'/input_files/dimensions.csv')
    
    
    ############################################################
    # save actual ratio_df values
    ############################################################
    collection = []
    index = []
    for f in files_use:
        gene_name = f.split('.')[0].split('_')[0]
        use_zerod_A = flags_df.loc[gene_name]['use_zerod_A_matrix']
        df_name = gene_name+'_zerod'+str(use_zerod_A)+'_mRNA_ratios_and_MA_vals.csv'
        ratios_df = pd.read_csv('../data/processed/saved_mRNA_ratios_MA_vals/'+df_name, index_col = 0)
        collection.append(ratios_df['actual_mRNA_ratio'])
        index.append(gene_name)
    ratios_combo_df = pd.DataFrame(collection, index = index)
    ratios_combo_df.T.to_csv(GAMS_run_dir+'/input_files/actual_mRNA_ratio.csv')
    

    ############################################################
    # create constants files for mRNA calculation
    ############################################################
    index = []
    collection = []
    for f in files_use:
        gene_name = f.split('.')[0].split('_')[0]
        gene_grid_name = '../data/interim/gene_grid_constants/'+gene_name+'.pkl'
        pickle_in = open(gene_grid_name, 'rb')
        grid_constants = pickle.load(pickle_in)
        pickle_in.close()
        index.append(gene_name)
        collection.append(grid_constants)
    constants_df = pd.DataFrame(collection, index = index)
    constants_df.T.to_csv(GAMS_run_dir+'/input_files/grid_constants.csv')     
    # make sample_constants file, for now just contains RNAP conc
    RNAP_conc_df = pd.read_csv('../data/interim/sample_constants/RNAP_conc.csv', index_col = 0)
    RNAP_conc_df = RNAP_conc_df.loc[ratios_df.index]
    RNAP_conc_df.T.to_csv(GAMS_run_dir+'/input_files/sample_constants.csv')     
    
    
    ############################################################
    # create TF concentration
    ############################################################
    concentrations = pd.read_csv('../data/external/validation_data_sets/converted_log_tpm_in_M.csv', index_col = 0)
    keep = []
    old_to_new = {}
    for promoter in set(flags_df.act_iM):
        if str(promoter) == 'nan':
            continue
        gene = iM_to_b_regulator[promoter][0] # not sure how we're handling this long term, for now there shall only be one regulator of each iModulon
        keep.append(gene)
        old_to_new.update({gene : promoter})
    concentrations.loc[keep].rename(index = old_to_new).to_csv(GAMS_run_dir+'/input_files/exported_act_TF_conc.csv')
    keep = []
    old_to_new = {}
    for inhibitor in set(flags_df.inh_iM):
        if str(inhibitor) == 'nan':
            continue
        gene = iM_to_b_regulator[inhibitor][0]
        keep.append(gene)
        old_to_new.update({gene : inhibitor})
    concentrations.loc[keep].rename(index = old_to_new).to_csv(GAMS_run_dir+'/input_files/exported_inh_TF_conc.csv')
    
    
    ############################################################
    # run GAMS
    ############################################################
    # call GAMS
    shutil.copyfile('../GAMS/conopt.opt', GAMS_run_dir+'/conopt.opt')
    if platform.system() == 'Windows':
        #gams_loc =  r'"C:\GAMS\win64\28.2\gams.exe"'
        shutil.copyfile('../GAMS/merged_iMs_model.gms', GAMS_run_dir+'/combined_model.gms')
        _ = subprocess.call(GAMS_exec+' combined_model.gms', cwd = GAMS_run_dir, shell=True)
    else:
        shutil.copyfile('../GAMS/merged_iMs_model.gms', GAMS_run_dir+'/combined_model.gms')
        if stable_flags['supress_output']:
            _ = subprocess.call(GAMS_exec+' combined_model > /dev/null', shell = True, cwd = GAMS_run_dir)
        else:
            _ = subprocess.call(GAMS_exec+' combined_model', shell = True, cwd = GAMS_run_dir)

def read_multi_GAMs(GAMS_run_dir):
    """
    Reads GAMS results from specific folder
    
    Inputs:
        GAMS_run_dir (string) : output location from a past GAMS run
    
    Returns:
        mRNA_ratio_df (dataframe) : mRNA ratio dataframe used as input for GAMS
        GAMS_cAct (dataframe) : cActivator values output from GAMS
        act_kd_df (dataframe) : activator Kd values output from GAMS
        act_metab (dataframe) : activator metabolite concentrations output from GAMS
        act_kd_metab_df (dataframe) : activator metabolite Kd values output from GAMS
        GAMS_cInh (dataframe) : cInhibitor values output from GAMS
        inh_kd_df (dataframe) : inhibitor Kd values output from GAMS
        inh_metab (dataframe) : inhibitor metabolite concentrations output from GAMS
        inh_kd_metab_df (dataframe) : inhibitor metabolite Kd values output from GAMS
        
    """
    # look at GAMs results
    
    # inputs
    saved_cActivators = pd.read_csv(GAMS_run_dir+'/input_files/composite_cAct_vals.csv', index_col = 0)
    saved_cInhibitors = pd.read_csv(GAMS_run_dir+'/input_files/composite_cInh_vals.csv', index_col = 0)
    act_TF_concs = pd.read_csv(GAMS_run_dir+'/input_files/exported_act_TF_conc.csv', index_col = 0)
    inh_TF_concs = pd.read_csv(GAMS_run_dir+'/input_files/exported_inh_TF_conc.csv', index_col = 0)
    input_constants = pd.read_csv(GAMS_run_dir+'/input_files/TF_constants.csv', index_col = 0)
    cAct_mapping = pd.read_csv(GAMS_run_dir+'/input_files/cAct_mapping.csv', index_col = 0)
    cInh_mapping = pd.read_csv(GAMS_run_dir+'/input_files/cInh_mapping.csv', index_col = 0)
    
    # outputs
    act_kd_df = pd.read_csv(GAMS_run_dir+'/output_files/cAct_Kd_results.csv', index_col = 0)
    inh_kd_df = pd.read_csv(GAMS_run_dir+'/output_files/cInh_Kd_results.csv', index_col = 0)
    act_metab = pd.read_csv(GAMS_run_dir+'/output_files/act_metab_Total.csv', index_col = 0)
    inh_metab = pd.read_csv(GAMS_run_dir+'/output_files/inh_metab_Total.csv', index_col = 0)

    # limit overlaps
    saved_cActivators = saved_cActivators.loc[act_kd_df.index]
    saved_cInhibitors = saved_cInhibitors.loc[inh_kd_df.index]

    
    # let's unlog kd_df and metab_df
    act_kd_df[act_kd_df.columns[1:]] = 10**act_kd_df[act_kd_df.columns[1:]].astype(float)
    inh_kd_df[inh_kd_df.columns[1:]] = 10**inh_kd_df[inh_kd_df.columns[1:]].astype(float)
    act_metab[act_metab.columns[1:]] = 10**act_metab[act_metab.columns[1:]].astype(float)
    inh_metab[inh_metab.columns[1:]] = 10**inh_metab[inh_metab.columns[1:]].astype(float)
    
    
    # activator read in
    temp = pd.read_csv(GAMS_run_dir+'/output_files/GAMS_cAct.csv', index_col = 0)
    GAMS_cAct = pd.DataFrame(index = list(set(temp.index)), columns = list(set(temp.gene)))
    for col in GAMS_cAct.columns:
        GAMS_cAct[col] = temp[temp.gene == col].loc[GAMS_cAct.index]['Val']
    
    # inhibitor read in
    temp = pd.read_csv(GAMS_run_dir+'/output_files/GAMS_cInh.csv', index_col = 0)
    GAMS_cInh = pd.DataFrame(index = list(set(temp.index)), columns = list(set(temp.gene)))
    for col in GAMS_cInh.columns:
        GAMS_cInh[col] = temp[temp.gene == col].loc[GAMS_cInh.index]['Val']
            
    # read in activator and inhibitor metabolite Kd's
    act_kd_metab_df = pd.read_csv(GAMS_run_dir+'/output_files/act_Kd_metab.csv', index_col = 0)
    inh_kd_metab_df = pd.read_csv(GAMS_run_dir+'/output_files/inh_Kd_metab.csv', index_col = 0)
    act_kd_metab_df[act_kd_metab_df.columns[1:]] = 10**act_kd_metab_df[act_kd_metab_df.columns[1:]].astype(float)
    inh_kd_metab_df[inh_kd_metab_df.columns[1:]] = 10**inh_kd_metab_df[inh_kd_metab_df.columns[1:]].astype(float)
    
    # read in mRNA, convert to log tpm, output both
    input_df = pd.read_csv(GAMS_run_dir+'/output_files/calculated_mRNA.csv', index_col = 0)
    mRNA_ratio_df = pd.DataFrame(index = list(set(input_df.index)), columns = list(set(input_df.gene)))
    for index, row in input_df.iterrows():
        mRNA_ratio_df.at[index, row['gene']] = row['Val']
    
    # return
    return(mRNA_ratio_df, GAMS_cAct, act_kd_df, act_metab, act_kd_metab_df, GAMS_cInh, inh_kd_df, inh_metab, inh_kd_metab_df)