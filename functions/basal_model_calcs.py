"""
Contains functions pertaining to the calculation of basal constants for the transcriptomic model.

Functions:
basal_values - Inputs an equation string, settings flags, and number of steps, outputs various possible basal solutions

"""

# imports
import sys
sys.path.insert(0, '../functions/')
from sympy import *
import pandas as pd
import promoter_solving_core as ps
import numpy as np

# In the future, when we move towards a promoter based model, this code will need to be replaced with promoter specific calculations
def basal_values(eq_str, flags, num_steps = 3):
    """
    Inputs an equation string, settings flags, and number of steps, outputs various possible basal solutions
    
    Inputs:
        eq_str (string) : sypmy-friendly equation
        flags (dict) : dictionary of settings flags and constants values
        num_steps (int) : number of iterations of solutions for the basal values, total number of solutions will be equal to num_steps**2
    
    Returns:
        grid_constants (dict) : dictionary of returned basal constants
    """

    # Define constants
    log_test = {
        'KdRNAP': [-8, -3], # [-7, -5]
        'kEscape': [-3, 1], # [-3, 1]
    }
    
    grid_constants = {
        #'KeqOpening': 10**-0.34444956947383365, gets set later
        'RNAP': flags['cell_constants']['RNAP'],#10**-6,
        'kDeg': flags['cell_constants']['kDeg'],#np.log(2)/t_half_life_deg, # Rate of degradation
        'promoterConcVal': flags['cell_constants']['promoterConcVal'],#10**-9, # Promoter concentration
        'TF': 0,#0, # Concentration of TF
        'u': flags['cell_constants']['u'],#1/3600, # Growth Rate
        'mRNA_total': flags['cell_constants']['mRNA_total'],
        'cell_volume': flags['cell_constants']['cell_volume'],
        'k_d_TF': 1, # May change depending on model
    }
    
    # Parameter Equation
    parameter_equation = sympify('Eq((KeqOpening*kEscape*promoterConcVal)/((KdRNAP/RNAP+1+KeqOpening+KdRNAP/RNAP*TF/k_d_TF)*(u+kDeg)),mRNA)')

    # Load in the precise data for gene expression
     # loading
    if flags['include_Amy_samples']:
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
    precise_data = log_tpm_df
        
    gene_exp = [2**precise_data.loc[flags['central_gene'], flags['basal_conditions']].mean(axis = 0)]

    # Only using num_steps = 3 to make it easier to manually iterate through the values, feel free to increase it if you'd like
    lambda_df, k_df = ps.create_grid(gene_exp = gene_exp, gene_name = [flags['central_gene']], equation = parameter_equation, constants = grid_constants, num_steps = num_steps, **log_test)
    
    grid_use = flags['grid_use']
    if grid_use == 'median':
        grid_use = int(len(k_df.iloc[0, 1]) / 2)
    grid_use = int(grid_use)
    
    grid_vals = k_df.iloc[0, 1][grid_use]
    grid_constants['KdRNAP'] = 10**(grid_vals[0])
    grid_constants['kEscape'] = 10**(grid_vals[1])
    grid_constants['KeqOpening'] = 10**(grid_vals[2])

    return(grid_constants)