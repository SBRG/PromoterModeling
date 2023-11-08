# basal model calculations
import sys
sys.path.insert(0, '../functions/')
from sympy import *
import pandas as pd
import promoter_solving_core as ps
import numpy as np


# In the future, when we move towards a promoter based model, this code will need to be replaced with promoter specific calculations
def basal_values(eq_str, flags):

    # Define constants
    log_test = {
        'KdRNAP': [-7,-5],
        'kEscape': [-3,1],
    }
    
    t_half_life_deg = 300
    grid_constants = {
        'KdRNAP': 10**-5,
        'KdRNAPCrp': 2.5118864315095796e-07*1.4,
        #'KeqOpening': 10**-0.34444956947383365, gets set later
        'RNAP': 10**-6,
        'mRNA_total': 1800, # Total mRNA/cell from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3554401
        'cell_volume': 10**-15, # Liters from https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100004&ver=19
        'k_d_TF': 1, # May change depending on model
        'kDeg': np.log(2)/t_half_life_deg, # Rate of degradation
        'promoterConcVal': 10**-9, # Promoter concentration
        'TF': 0, # Concentration of TF
        'u': 1/3600, # Growth Rate
    }
    
    # Parameter Equation
    parameter_equation = sympify('Eq((KeqOpening*kEscape*promoterConcVal)/((KdRNAP/RNAP+1+KeqOpening+KdRNAP/RNAP*TF/k_d_TF)*(u+kDeg)),mRNA)')

    # Load in the precise data for gene expression
    # NOTE: In the basal model I'm building, I am using the Precise1k data, but in this cell I will try to use the Precise1.0
    precise_path = '../data/precise_1.0/log_tpm.csv'
    precise_data = pd.read_csv(filepath_or_buffer = precise_path, index_col = 'Unnamed: 0')

    gene_exp = [2**precise_data.loc[flags['central_gene'], flags['basal_conditions']].mean(axis = 0)]

    # Only using num_steps = 3 to make it easier to manually iterate through the values, feel free to increase it if you'd like
    lambda_df, k_df = ps.create_grid(gene_exp = gene_exp, gene_name = [flags['central_gene']], equation = parameter_equation, constants = grid_constants, num_steps = 3, **log_test)
    grid_vals = k_df.iloc[0, 1][-1]
    grid_constants['KdRNAP'] = 10**(grid_vals[0])
    grid_constants['kEscape'] = 10**(grid_vals[1])
    grid_constants['KeqOpening'] = 10**(grid_vals[2])

    return(grid_constants)