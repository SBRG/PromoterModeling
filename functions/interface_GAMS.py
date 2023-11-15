# run GAMs, need to test, move to outside .py file
import subprocess
import os
import pandas as pd
import numpy as np
import pickle
import shutil

def run_GAMs(flags_df, stable_flags, promoter, inhibitor, cell_constants, GAMs_run_dir, parameter_flags = None):
    ############################################################
    # save input parameters
    ############################################################
    if not parameter_flags:
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
    }
    df = pd.DataFrame(list(baby_dict.items()), columns=['Parameter', 'Value']).set_index('Parameter')
    df.to_csv(GAMs_run_dir+'/input_files/parameters.csv')
    
    
    ############################################################
    # create TF concentration
    ############################################################
    if stable_flags['include_Amy_samples']:
        # merge together log_tpm_df files
        log_tpm_df = pd.read_csv('../data/precise_1.0/log_tpm.csv', index_col = 0)
        starve_log_tpm = pd.read_csv('../data/validation_data_sets/stationary_phase/cleaned_log_tpm_qc.csv', index_col = 0)
        to_blank_inds = list(set(log_tpm_df.index) - set(starve_log_tpm.index))
        # need to create zero rows for missing values
        zeros_data = {col : 0 for col in starve_log_tpm.columns}
        zeros_df = pd.DataFrame(zeros_data, index = to_blank_inds)
        starve_log_tpm = pd.concat([starve_log_tpm, zeros_df])
        starve_log_tpm = starve_log_tpm.loc[log_tpm_df.index]
        log_tpm_df = pd.concat([starve_log_tpm, log_tpm_df], axis = 1)
    else:
        log_tpm_df = pd.read_csv('../data/precise_1.0/log_tpm.csv', index_col = 0)
    concentrations = 2**log_tpm_df*cell_constants['mRNA_total']/cell_constants['cell_volume']/(6.022*10**23)
    concentrations.loc[promoter].to_csv(GAMs_run_dir+'/input_files/exported_act_TF_conc.csv')
    concentrations.loc[inhibitor].to_csv(GAMs_run_dir+'/input_files/exported_inh_TF_conc.csv')

    
    # first let's merge together the files
    files_use = []
    if stable_flags['run_on_all']:
        files = os.listdir('../data/cAct_cInh_vals/')
        for f in files:
            gene_name = file.split('.')[0].split('_')[0]
            if stable_flags['use_greedy'] and 'greed' not in f: continue
            if not stable_flags['use_greedy'] and 'greed' in f: continue
            files_use.append(f)
    else:
        for sample in stable_flags['limit_samples']:
            use_zerod_A = flags_df.loc[sample]['use_zerod_A_matrix']
            if stable_flags['use_greedy']:
                f_name = sample+'_greedy.pkl'
            else:
                f_name = sample+'.pkl'
            if os.path.exists('../data/cAct_cInh_vals/'+f_name):
                files_use.append(f_name)
    if False: # old way
        if stable_flags['run_on_all']:
            files = os.listdir('../data/save_for_GAMs/')
            for f in files:
                gene_name = file.split('_')[0]
                use_zerod_A = flags_df.loc[gene_name]['use_zerod_A_matrix']
                if 'composite' in f: continue
                if 'zerod'+str(use_zerod_A) not in f: continue
                if stable_flags['use_greedy'] and 'non_greed' in f: continue
                if not stable_flags['use_greedy'] and 'non_greed' not in f: continue
                files_use.append(f)
        else:
            for sample in stable_flags['limit_samples']:
                use_zerod_A = flags_df.loc[sample]['use_zerod_A_matrix']
                if stable_flags['use_greedy']:
                    f_name = sample+'_zerod'+str(use_zerod_A)+'_cAct_cInh_vals.csv'
                else:
                    f_name = sample+'_non_greed_zerod'+str(use_zerod_A)+'_cAct_cInh_vals.csv'
                if os.path.exists('../data/save_for_GAMs/'+f_name):
                    files_use.append(f_name)
             
                
    for f in files_use:
        shared_indices = []
        first = True
        for file in files_use:
            gene_name = file.split('.')[0].split('_')[0]
            temp_df = pd.read_pickle('../data/cAct_cInh_vals/'+file)
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
            act_df[gene_name] = pd.read_pickle('../data/cAct_cInh_vals/'+file)['cAct'].loc[act_df.index].values
            inh_df[gene_name] = pd.read_pickle('../data/cAct_cInh_vals/'+file)['cInh'].loc[inh_df.index].values
        
    act_df.to_csv(GAMs_run_dir+'/input_files/composite_cAct_vals.csv')
    inh_df.to_csv(GAMs_run_dir+'/input_files/composite_cInh_vals.csv')
    vals = []
    for col in act_df.columns:
        vals.append(max(act_df[col]))
    max_df = pd.DataFrame(vals, index = act_df.columns)
    max_df.to_csv(GAMs_run_dir+'/input_files/max_cActs.csv')
    vals = []
    for col in inh_df.columns:
        vals.append(max(inh_df[col]))
    max_df = pd.DataFrame(vals, index = inh_df.columns)
    max_df.to_csv(GAMs_run_dir+'/input_files/max_cInhs.csv')
    
    ############################################################
    # create constants file for mRNA calculation
    ############################################################
    index = []
    collection = []
    for f in files_use:
        gene_name = f.split('.')[0].split('_')[0]
        gene_grid_name = '../data/gene_grid_constants/'+gene_name+'.pkl'
        pickle_in = open(gene_grid_name, 'rb')
        grid_constants = pickle.load(pickle_in)
        pickle_in.close()
        index.append(gene_name)
        collection.append(grid_constants)
    constants_df = pd.DataFrame(collection, index = index)
    constants_df.T.to_csv(GAMs_run_dir+'/input_files/grid_constants.csv')     
    
    ############################################################
    # save actual ratio_df values
    ############################################################
    collection = []
    index = []
    for f in files_use:
        gene_name = f.split('.')[0].split('_')[0]
        use_zerod_A = flags_df.loc[gene_name]['use_zerod_A_matrix']
        df_name = gene_name+'_zerod'+str(use_zerod_A)+'_mRNA_ratios_and_MA_vals.csv'
        ratios_df = pd.read_csv('../data/saved_mRNA_ratios_MA_vals/'+df_name, index_col = 0)
        collection.append(ratios_df['actual_mRNA_ratio'])
        index.append(gene_name)
    ratios_combo_df = pd.DataFrame(collection, index = index)
    ratios_combo_df.T.to_csv(GAMs_run_dir+'/input_files/actual_mRNA_ratio.csv')
    
    # remove old results
    if stable_flags['delete_old']:
        if os.path.exists('../data/GAMS_output/cInh_Kd_results.csv'):
            os.remove('../data/GAMS_output/cInh_Kd_results.csv')
        if os.path.exists('../data/GAMS_output/cInh_TF_conc_results.csv'):
            os.remove('../data/GAMS_output/cInh_TF_conc_results.csv')
        if os.path.exists('../data/GAMS_output/cAct_Kd_results.csv'):
            os.remove('../data/GAMS_output/cAct_Kd_results.csv')
        if os.path.exists('../data/GAMS_output/cAct_TF_conc_results.csv'):
            os.remove('../data/GAMS_output/cAct_TF_conc_results.csv')

    # call GAMs
    if stable_flags['run_seperate']:
        shutil.copyfile('../GAMs/cAct_model.gms', GAMs_run_dir+'/cAct_model.gms')
        shutil.copyfile('../GAMs/cInh_model.gms', GAMs_run_dir+'/cInh_model.gms')
        if stable_flags['supress_output']:
            _ = subprocess.call('gams cAct_model > /dev/null', shell = True, cwd = GAMs_run_dir)
            _ = subprocess.call('gams cInh_model > /dev/null', shell = True, cwd = GAMs_run_dir)
        else:
            _ = subprocess.call('gams cAct_model', shell = True, cwd = GAMs_run_dir)
            _ = subprocess.call('gams cInh_model', shell = True, cwd = GAMs_run_dir)
    else:
        shutil.copyfile('../GAMs/combined_model.gms', GAMs_run_dir+'/combined_model.gms')
        if stable_flags['supress_output']:
            _ = subprocess.call('gams combined_model > /dev/null', shell = True, cwd = GAMs_run_dir)
        else:
            _ = subprocess.call('gams combined_model', shell = True, cwd = GAMs_run_dir)
    
    
    
    
    
    
def read_GAMs(GAMs_run_dir):
    # look at GAMs results

    # load in cActivators
    saved_cActivators = pd.read_csv(GAMs_run_dir+'/input_files/composite_cAct_vals.csv', index_col = 0)

    # GAMS calculated cActivators
    cAct_kd_df = 10**pd.read_csv(GAMs_run_dir+'/output_files/cAct_Kd_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators[cAct_kd_df.columns]
    cAct_TF_conc_df = 10**pd.read_csv(GAMs_run_dir+'/output_files/cAct_TF_conc_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators.loc[cAct_TF_conc_df.columns]
    calc_cAct = pd.DataFrame(index = saved_cActivators.columns, columns = saved_cActivators.index)
    cActs = []
    for sample in calc_cAct.columns:
        for gene in calc_cAct.index:
            calc_cAct.at[gene, sample] = cAct_TF_conc_df[sample].values[0] / cAct_kd_df[gene].values[0]
    calc_cAct = calc_cAct.T


    # now cInhibitor
    # load in cActivators
    saved_cActivators = pd.read_csv(GAMs_run_dir+'/input_files/composite_cInh_vals.csv', index_col = 0)

    # GAMS calculated cActivators
    kd_df = 10**pd.read_csv(GAMs_run_dir+'/output_files/cInh_Kd_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators[kd_df.columns]
    TF_conc_df = 10**pd.read_csv(GAMs_run_dir+'/output_files/cInh_TF_conc_results.csv', index_col = 0).astype(float).T
    saved_cActivators = saved_cActivators.loc[TF_conc_df.columns]
    calc_cInh = pd.DataFrame(index = saved_cActivators.columns, columns = saved_cActivators.index)
    cActs = []
    for sample in calc_cInh.columns:
        for gene in calc_cInh.index:
            calc_cInh.at[gene, sample] = TF_conc_df[sample].values[0] / kd_df[gene].values[0]
    calc_cInh = calc_cInh.T
    
    # return
    return(calc_cAct, cAct_kd_df, cAct_TF_conc_df, calc_cInh, kd_df, TF_conc_df)