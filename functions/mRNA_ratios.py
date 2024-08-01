"""
Creates mRNA ratio values and MA values for a specific gene

Functions:
calculate_mRNA_ratios_and_MA_values - Inputs iModulons and flags to create mRNA ratio values and MA values for a specific gene

"""

# imports
import pandas as pd
import numpy as np
import pickle
import os

def calculate_mRNA_ratios_and_MA_values(iM_act, iM_inh, input_parameters):
    """
    Inputs an equation string, settings flags, and number of steps, outputs various possible basal solutions
    
    Inputs:
        iM_act (string) : activating iModulon
        iM_inh (string) : inhibiting iModulon
        input_parameters (int) : dictionary of settings flags and constants values, contains central gene to run on
    
    Returns:
        values_df (dataframe) : contains mRNA ratios and MA values for a specific gene
    """
        
    # unload flags
    gene = input_parameters['central_gene']
    use_zerod_A_matrix = input_parameters['use_zerod_A_matrix']
    basal_conditions = input_parameters['basal_conditions']
    
    
    # loading
    if input_parameters['include_Amy_samples']:
        # merge together log_tpm_df files
        log_tpm_df = pd.read_csv('../data/external/imodulon_info/log_tpm.csv', index_col = 0)
        starve_log_tpm = pd.read_csv('../data/external/validation_data_sets/stationary_phase/cleaned_log_tpm_qc.csv', index_col = 0)
        to_blank_inds = list(set(log_tpm_df.index) - set(starve_log_tpm.index))
        # need to create zero rows for missing values
        zeros_data = {col : 0 for col in starve_log_tpm.columns}
        zeros_df = pd.DataFrame(zeros_data, index = to_blank_inds)
        starve_log_tpm = pd.concat([starve_log_tpm, zeros_df])
        starve_log_tpm = starve_log_tpm.loc[log_tpm_df.index]
        log_tpm_df = pd.concat([starve_log_tpm, log_tpm_df], axis = 1)
    else:
        log_tpm_df = pd.read_csv('../data/external/imodulon_info/log_tpm.csv', index_col = 0)
    
    iM_table = pd.read_csv('../data/external/imodulon_info/iM_table.csv', index_col = 0)
    k_to_iM = {index : name for index, name in zip(iM_table.index, iM_table.name)}
    M_df = pd.read_csv('../data/external/imodulon_info/M.csv', index_col = 0).rename(columns = {str(k) : v for k, v in k_to_iM.items()})
    
    overlap = list(set(M_df.index).intersection(log_tpm_df.index))
    M_df = M_df.loc[overlap]
    log_tpm_df = log_tpm_df.loc[overlap]
        

    #iM_table = pd.read_csv('../data/precise_1.0/iM_table.csv', index_col = 0)
    #M_df = M_df.rename(columns = {str(index) : row['name'] for index, row in iM_table.iterrows()})
    
    # filter out samples set by paramater_optimization/6_find_outliers_to_drop.ipynb
    if input_parameters['remove_outliers']:
        dict_path = '../data/interim/misc_dictionaries/case_to_mRNA_passed.pkl'
        if os.path.exists(dict_path):
            pickle_in = open(dict_path, 'rb')
            case_to_mRNA_passed = pickle.load(pickle_in)
            pickle_in.close()
        else:
            case_to_mRNA_passed = {}
        
        case = '__'.join([str(iM) for iM in [iM_act, iM_inh]]).replace(' ', '_').replace('/', '_')
        if case in case_to_mRNA_passed:
            passed = case_to_mRNA_passed[case]
            for cond in basal_conditions:
                passed.append(cond)
            log_tpm_df = log_tpm_df[list(set(passed).intersection(set(log_tpm_df.columns)))]
    
    # creates zerod matrices
    if use_zerod_A_matrix:
        gene_iMs_df = pd.read_csv('../data/external/imodulon_info/corrected/gene_presence_matrix.csv', index_col = 0)
        gene_iMs_df.columns = M_df.columns
        if type(iM_act) == float:
            genes_to_zero = list(gene_iMs_df.index[[val for val in gene_iMs_df[[iM_inh]].T.any()]])
            iMs_to_zero = list(set(gene_iMs_df.columns) - set([iM_inh]))
        elif type(iM_inh) == float:
            genes_to_zero = list(gene_iMs_df.index[[val for val in gene_iMs_df[[iM_act]].T.any()]])
            iMs_to_zero = list(set(gene_iMs_df.columns) - set([iM_act]))
        else:
            genes_to_zero = list(gene_iMs_df.index[[val for val in gene_iMs_df[[iM_act, iM_inh]].T.any()]])
            iMs_to_zero = list(set(gene_iMs_df.columns) - set([iM_act, iM_inh]))

        zerod_M = M_df.copy()
        zerod_M.loc[genes_to_zero, iMs_to_zero] = 0
        #zerod_M = zerod_M.drop(columns = ['fps__fps_ptsI_ale3__1', 'fps__fps_ptsI_ale3__2', 'fps__fps_ptsI_ale1__1', 'fps__fps_ptsI_ale1__2'])

        # Calculate the inverse of DataFrame M
        M_inverse = pd.DataFrame(np.linalg.pinv(zerod_M.values), zerod_M.columns, zerod_M.index)

        # Solve for DataFrame A: A = M_inverse * X
        #fixed_X = log_tpm_df.div(log_tpm_df[basal_conditions].mean(axis = 1), axis = 'index')
        if input_parameters['basal_or_hard_val'] == 'basal':
            fixed_X = log_tpm_df.sub(log_tpm_df[basal_conditions].mean(axis = 1), axis = 'index')
        else:
            fixed_X = log_tpm_df.sub(input_parameters['hard_val'], axis = 'index')
        to_drop = list(set(['fps__fps_ptsI_ale3__1', 'fps__fps_ptsI_ale3__2', 'fps__fps_ptsI_ale1__1', 'fps__fps_ptsI_ale1__2']).intersection(set(fixed_X.columns))) # these samples are problematic
        fixed_X = fixed_X.fillna(0).drop(columns = to_drop)
        zerod_A_df = M_inverse.dot(fixed_X)

        A_df = zerod_A_df
        M_df = zerod_M
    
    # create data matrix
    act_MAs = []
    inh_MAs = []
    ratio_act_MAs = []
    ratio_inh_MAs = []
    index = []
    actual_counts = []
    if input_parameters['basal_or_hard_val'] == 'basal':
        log_x_c = log_tpm_df.loc[gene][basal_conditions].mean()
    else:
        log_x_c = input_parameters['hard_val']
    p1k_center = log_tpm_df.loc[gene][['p1k_00001', 'p1k_00002']].mean()
    
    # predict mRNA values
    for key, _ in (A_df.loc[A_df.index[0]].T*(M_df[M_df.columns[0]].loc[gene])).items():
        index.append(key)
        actual_counts.append(2**(log_tpm_df.loc[gene][key]) / 2**(log_x_c))
    if type(iM_act) == float:
        act_MAs = [0 for _ in index]
        ratio_act_MAs = [0 for _ in index]
    else:
        for key, val in (A_df.loc[iM_act].T*(M_df[iM_act].loc[gene])).items():
            act_MAs.append(val)
            ratio_act_MAs.append(2**(val+p1k_center) / 2**log_x_c)
    if type(iM_inh) == float:
        inh_MAs = [0 for _ in index]
        ratio_inh_MAs = [0 for _ in index]
    else:
        for key, val in (A_df.loc[iM_inh].T*(M_df[iM_inh].loc[gene])).items():
            inh_MAs.append(val)
            ratio_inh_MAs.append(2**(val+p1k_center) / 2**log_x_c)

    values_df = pd.DataFrame(index = index)
    values_df['MA_activator'] = act_MAs
    values_df['MA_inhibitor'] = inh_MAs
    values_df['actual_mRNA_ratio'] = actual_counts
    values_df['ratio_MA_activator'] = ratio_act_MAs
    values_df['ratio_MA_inhibitor'] = ratio_inh_MAs
    
    # if set, remove basal conditions from sample to not weight results towards them
    if input_parameters['drop_basal_conds']:
        new_group = list(values_df.index)
        for cond in basal_conditions:
            new_group.remove(cond)
        values_df = values_df.loc[new_group]

        
    return(values_df)