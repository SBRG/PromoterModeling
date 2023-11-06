# run GAMs, need to test, move to outside .py file
import subprocess
import os
import pandas as pd
import numpy as np


def run_GAMs(flags, regulator, cell_constants):
    ############################################################
    # create TF concentration
    ############################################################
    log_tpm_df = pd.read_csv('../data/precise_1.0/log_tpm.csv', index_col = 0)
    concentrations = 2**log_tpm_df*10**-6*cell_constants['mRNA_total']/cell_constants['cell_volume']/(6.022*10**23)
    concentrations.loc[regulator].to_csv('../data/save_for_GAMs/exported_TF_conc.csv')
    
    
    # first let's merge together the files
    files_use = []
    if flags['run_on_all']:
        files = os.listdir('../data/save_for_GAMs/')
        for f in files:
            if 'composite' in f: continue
            if 'zerod'+str(flags['use_zerod_A_matrix']) not in f: continue
            if flags['use_greedy'] and 'non_greed' in f: continue
            if not flags['use_greedy'] and 'non_greed' not in f: continue
            files_use.append(f)
    else:
        for sample in flags['limit_samples']:
            if flags['use_greedy']:
                f_name = sample+'_zerod'+str(flags['use_zerod_A_matrix'])+'_cAct_cInh_vals.csv'
            else:
                f_name = sample+'_non_greed_zerod'+str(flags['use_zerod_A_matrix'])+'_cAct_cInh_vals.csv'
            if os.path.exists('../data/save_for_GAMs/'+f_name):
                files_use.append(f_name)
    for f in files_use:
        shared_indices = []
        first = True
        for file in files_use:
            gene_name = file.split('_')[0]
            temp_df = pd.read_csv('../data/save_for_GAMs/'+file, index_col = 0)
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
            gene_name = file.split('_')[0]
            act_df[gene_name] = pd.read_csv('../data/save_for_GAMs/'+file, index_col = 0)['cAct'].loc[act_df.index].values
            inh_df[gene_name] = pd.read_csv('../data/save_for_GAMs/'+file, index_col = 0)['cInh'].loc[inh_df.index].values
    act_df.to_csv('../data/save_for_GAMs/composite_cAct_vals.csv')
    inh_df.to_csv('../data/save_for_GAMs/composite_cInh_vals.csv')

    # remove old results
    if flags['delete_old']:
        if os.path.exists('../data/GAMS_output/cInh_Kd_results.csv'):
            os.remove('../data/GAMS_output/cInh_Kd_results.csv')
        if os.path.exists('../data/GAMS_output/cInh_TF_conc_results.csv'):
            os.remove('../data/GAMS_output/cInh_TF_conc_results.csv')
        if os.path.exists('../data/GAMS_output/cAct_Kd_results.csv'):
            os.remove('../data/GAMS_output/cAct_Kd_results.csv')
        if os.path.exists('../data/GAMS_output/cAct_TF_conc_results.csv'):
            os.remove('../data/GAMS_output/cAct_TF_conc_results.csv')

    # call GAMs
    _ = subprocess.call('gams cAct_model', shell = True, cwd = '../GAMs')
    _ = subprocess.call('gams cInh_model', shell = True, cwd = '../GAMs')
    
def read_GAMs(flags):
    # look at GAMs results

    # load in cActivators
    saved_cActivators = pd.read_csv('../data/save_for_GAMs/composite_cAct_vals.csv', index_col = 0)

    # GAMS calculated cActivators
    cAct_kd_df = 10**pd.read_csv('../data/GAMS_output/cAct_Kd_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators[cAct_kd_df.columns]
    cAct_TF_conc_df = 10**pd.read_csv('../data/GAMS_output/cAct_TF_conc_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators.loc[cAct_TF_conc_df.columns]
    calc_cAct = pd.DataFrame(index = saved_cActivators.columns, columns = saved_cActivators.index)
    cActs = []
    for sample in calc_cAct.columns:
        for gene in calc_cAct.index:
            calc_cAct.at[gene, sample] = cAct_TF_conc_df[sample].values[0] / cAct_kd_df[gene].values[0]
    calc_cAct = calc_cAct.T


    # now cInhibitor
    # load in cActivators
    saved_cActivators = pd.read_csv('../data/save_for_GAMs/composite_cInh_vals.csv', index_col = 0)

    # GAMS calculated cActivators
    kd_df = 10**pd.read_csv('../data/GAMS_output/cInh_Kd_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators[kd_df.columns]
    TF_conc_df = 10**pd.read_csv('../data/GAMS_output/cInh_TF_conc_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators.loc[TF_conc_df.columns]
    calc_cInh = pd.DataFrame(index = saved_cActivators.columns, columns = saved_cActivators.index)
    cActs = []
    for sample in calc_cInh.columns:
        for gene in calc_cInh.index:
            calc_cInh.at[gene, sample] = TF_conc_df[sample].values[0] / kd_df[gene].values[0]
    calc_cInh = calc_cInh.T
    
    # return
    return(calc_cAct, cAct_kd_df, cAct_TF_conc_df, calc_cInh, kd_df, TF_conc_df)