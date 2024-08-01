"""
Contains functions pertaining to the calculation of basal constants for the transcriptomic model.

Functions:
create_lambdas - Griff Hughes made this function and I don't feel qualified to write his descriptions for him.
evaluate_lambda - Griff Hughes made this function and I don't feel qualified to write his descriptions for him.
create_shared_lambda_df - Griff Hughes made this function and I don't feel qualified to write his descriptions for him.
mRNA_cInhibitor_to_cActivator - Converts mRNA ratio and cInhibitor to a cActivator value
mRNA_cActivator_to_cInhibitor - Converts mRNA ratio and cActivator to a cInhibitor value
cActivator_cInhibitor_to_mRNA - Converts cInhibitor and cActivator to a mRNA ratio value
pick_KdRNAPActivator - Selects an optimal KdRNAPAct value
"""

# imports
import pandas as pd
import pickle
import numpy as np
from scipy.optimize import minimize
from sympy import *
from itertools import product
import math
import matplotlib.pyplot as plt


# these values should get overwritten in the actual code, but need to be initialized
KdRNAP = None
KeqOpening = None
RNAP = None
KdRNAPAct = None
lambda_df = None

# this is taken from promoter_solving_core.py to avoid cyclical imports, they should technically just both import
# it from a third location but eh, something to deal with later
# Modules

# Functions
def create_lambdas(
        equation: str,
        constants: dict
    ) -> pd.DataFrame:
    """
    Creates a df where each variable has a corresponding lambda function that 
        can be used to calculate the variable
    
    :param str equation: String of the form 'Eq(x,y)' that sympify can parse
    :param dict constants: Dict where the keys are symbolic variables and the 
        value is the constant's value
    """

    # Create symbolic variables from the input equation
    eqn = sympify(equation)

    # Substitute in constants
    eqn = eqn.subs(constants.items())

    # Identify parameters and sort them since free_symbols returns a set
    parameters = sorted(list(str(parameter) for parameter in eqn.free_symbols),
                        key = str.lower)

    # Create df to store parameter equations, lambda functions, and tuples
    lambda_df = pd.DataFrame(columns = ['equation',
                                        'lambda',
                                        'order'],
                             index =  parameters)

    for parameter in parameters:
        # Algebraically solve the equation for the parameter and apply constants
        lambda_df.loc[parameter,'equation'] = solve(eqn, parameter)

        # Create a tuple of parameters to sub into the lambdified equation
        lambda_df.loc[parameter,'order'] = tuple([p for p in parameters if p != parameter])

        # Lambdify each parameter's equation
        lambda_df.loc[parameter,'lambda'] = lambdify([lambda_df.loc[parameter,'order']],
                                                     lambda_df.loc[parameter,'equation'])

    return lambda_df

def evaluate_lambda(
        solve: str,
        lambda_df: pd.DataFrame,
        values: dict
    ) -> float:
    """
    Evaluates the solution to the promoter equation

    :param str solve: Parameter to solve for
    :param pd.DataFrame lambda_df: df containing the lambda functions and the 
        parameter order for the tuple to pass into the lambda function
    :param dict values: Key is parameter, val is the log10 of the value
    """

    # Create a tuple in the correct order to pass into the lambda function
    values_tuple = tuple([values[p] for p in lambda_df.loc[solve,'order']])

    # Evaluate the lambda function
    return((lambda_df.loc[solve,'lambda'](values_tuple))[0])

def create_shared_lambda_df(equation_string, grid):
    global lambda_df, KdRNAP, KeqOpening, RNAP
    our_grid = grid.copy()
    equation = sympify(equation_string)
    if 'KdRNAPAct' in our_grid:
        del(our_grid['KdRNAPAct'])
        
    # Create lambda functions that we can plug in to
    lambda_df = create_lambdas(equation, our_grid)
    
    KdRNAP = our_grid['KdRNAP']
    KeqOpening = our_grid['KeqOpening']
    RNAP = our_grid['RNAP']

def mRNA_cInhibitor_to_cActivator(mRNA, cInhibitor, KdRNAPAct, lambda_df_input = None):
    """
    Inputs mRNA ratio, cInhibitor and KdRNAPAct to predict cActivator
    
    Inputs:
        mRNA (float) : mRNA ratio value
        cInhibitor (float) : cInhibitor value
        KdRNAPAct (float) : KdRNAPAct value
        lambda_df_input (dataframe) : sympy lambda dataframe used for calculations
    
    Returns:
        cActivator (float) : mathematically derived cActivator value based on inputs
    """
    
    if type(lambda_df_input) == type(None):
        lambda_df_input = lambda_df
    cActivator = evaluate_lambda('cActivator', lambda_df_input, {'cInhibitor': cInhibitor, 'mRNARatio' : mRNA, 'KdRNAPAct' : KdRNAPAct})
    return(cActivator)

def mRNA_cActivator_to_cInhibitor(mRNA, cActivator, KdRNAPAct, lambda_df_input = None):
    """
    Inputs mRNA ratio, cActivator and KdRNAPAct to predict cInhibitor
    
    Inputs:
        mRNA (float) : mRNA ratio value
        cActivator (float) : cActivator value
        KdRNAPAct (float) : KdRNAPAct value
        lambda_df_input (dataframe) : sympy lambda dataframe used for calculations
    
    Returns:
        cInhibitor (float) : mathematically derived cInhibitor value based on inputs
    """
        
    if type(lambda_df_input) == type(None):
        lambda_df_input = lambda_df
    cInhibitor = evaluate_lambda('cInhibitor', lambda_df_input, {'cActivator': cActivator, 'mRNARatio' : mRNA, 'KdRNAPAct' : KdRNAPAct})
    return(cInhibitor)

def cActivator_cInhibitor_to_mRNA(cActivator, cInhibitor, KdRNAPAct, lambda_df_input = None):
    """
    Inputs cInhibitor, cActivator and KdRNAPAct to predict cInhibitor
    
    Inputs:
        cActivator (float) : cActivator value
        cInhibitor (float) : cInhibitor value
        KdRNAPAct (float) : KdRNAPAct value
        lambda_df_input (dataframe) : sympy lambda dataframe used for calculations
    
    Returns:
        mRNA (float) : mathematically derived mRNA ratio value based on inputs
    """
    
    if type(lambda_df_input) == type(None):
        lambda_df_input = lambda_df
    mRNA = evaluate_lambda('mRNARatio', lambda_df_input, {'cActivator': cActivator, 'cInhibitor' : cInhibitor, 'KdRNAPAct' : KdRNAPAct})
    return(mRNA)

def pick_KdRNAPActivator(ratios_df, flags):
    """
    Selects an optimal KdRNAPActivator value which optimizes two criteria in succession:
    1 - Finds the highest KdRNAPActivator value below KdRNAP which creates non-negative cActivator values
    2 - Selects a KdRNAPActivator value below said highest possible value which maximizes the spread of values to prevent the creation of extremely high cActivator outliers.
    
    Inputs:
        ratios_df (dataframe) : mRNA ratios dataframe
        flags (dict) : dictionary of settings flags and constants values
    
    Returns:
        optimal_KdRNAPAct (float) : optimal KdRNAPAct value
        ret_figs (array) : set of sanity plots used to verify
    """
    
    ret_figs = []
    initial_guess_ratio = flags['initial_guess_ratio']
    min_cInh = flags['base_cInhibitor_val']
    # v1 assumed cInhibitor = 0, this instead sets it to a minimum value instead of zero, I'm hoping this elevates values that are below zero and is a more realistic assumption
    
    # first let's find the maximum KdRNAPAct value for the mRNA values
    # using ratios_df is not good sampling, moving to a set spacing
    ratio_values = np.linspace(min(ratios_df['actual_mRNA_ratio']), max(ratios_df['actual_mRNA_ratio']), 1000)
    one_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - 1))
    def objective_function1(KdRNAPAct_temp):
        
        cActs = [mRNA_cInhibitor_to_cActivator(rat_val, min_cInh, KdRNAPAct_temp) for rat_val in ratio_values]

        # Set up penalty
        penalty = 0

        # I want to find the sample with the latest point above mRNA ratio 1 that has a negative cActivator value
        # i.e. I'm trying to figure out the highest KdRNAPAct value I can use that doesn't require negatives

        # let's search for the first index that has negative values
        ind = one_index + 1
        found = False
        while ind < len(cActs):
            if cActs[ind] < 0:
                # I've found it, I want to penalize it being very early, so I'll subtract how far along it is
                penalty -= (ind - one_index)
                found = True
                break
            ind += 1
        if not found: # to slope the left part towards the minimum
            penalty -= 1000000*KdRNAPAct_temp
        else:
            pass
        return(penalty)

    # setup and minimize
    initial_guess = KdRNAP * initial_guess_ratio
    bounds = [(1e-9, KdRNAP/.99)]
    
    result1 = minimize(objective_function1, initial_guess, method = 'Nelder-Mead', bounds = bounds)#, tol = 1e-10)
    max_KdRNAPAct = result1.x[0] 
    if flags['KdRNAPAct_sanity']:
        # Define the range of KdRNAPAct values for plotting
        KdRNAPAct_values = np.linspace(1e-9, KdRNAP / .99, 1000)  # Adjust the range as needed

        # Calculate the objective function values for each KdRNAPAct value
        objective_values = [objective_function1(KdRNAPAct) for KdRNAPAct in KdRNAPAct_values]

        # Create a plot to visualize the objective function
        fig, axs = plt.subplots(1, 2, figsize = (6, 3))
        axs[0].plot(KdRNAPAct_values, objective_values)
        axs[0].axvline(x = max_KdRNAPAct, c = 'k', ls = '--')
        #plt.ylim(min(objective_values), 10)
        axs[0].set_xlabel('KdRNAPAct')
        axs[0].set_ylabel('Objective Function Value')
        axs[0].set_title('1st Objective Function vs. KdRNAPAct')
    
    
    # now with our new maximum, let's look for the ideal value (and we can ignore negative value penalty I think)
    # Define your objective function
    max_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - max(ratio_values)))
    closest_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - .70*max(ratio_values)))
    def objective_function2(KdRNAPAct_temp):
        KdRNAPAct_temp = float(KdRNAPAct_temp)
        cActs = [mRNA_cInhibitor_to_cActivator(rat_val, min_cInh, KdRNAPAct_temp) for rat_val in ratio_values]

        # Set up penalty
        penalty = 0

        # let's try to maximize the 80% value
        penalty -= cActs[closest_index]

        # I want to maximize "spread" and getting a more diverse range of values... let me think on what that measure would be
        penalty -= np.corrcoef(cActs, np.arange(len(cActs)))[0, 1]
        penalty -= cActs[closest_index]
        
        if flags['use_target_range']:
            penalty += abs(10**flags['target_range'][1] - max(cActs))
            penalty += abs(10**flags['target_range'][0] - min(cActs))

        return(penalty)

    # setup and minimize
    initial_guess = max_KdRNAPAct / 2
    bounds = [(max_KdRNAPAct / 100, max_KdRNAPAct)]
    result2 = minimize(objective_function2, initial_guess, method = 'Nelder-Mead', bounds = bounds)#, tol = 1e-10)

    # The optimal KdRNAPAct value
    optimal_KdRNAPAct = result2.x[0]
    
    if flags['KdRNAPAct_sanity']:
        # Define the range of KdRNAPAct values for plotting
        KdRNAPAct_values = np.linspace(max_KdRNAPAct / 100, max_KdRNAPAct, 1000)  # Adjust the range as needed

        # Calculate the objective function values for each KdRNAPAct value
        objective_values = [objective_function2(KdRNAPAct) for KdRNAPAct in KdRNAPAct_values]

        # Create a plot to visualize the objective function  
        axs[1].plot(KdRNAPAct_values, objective_values)
        axs[1].axvline(x = optimal_KdRNAPAct, c = 'k', ls = '--')
        #plt.ylim(min(objective_values), 10)
        axs[1].set_xlabel('KdRNAPAct')
        axs[1].set_ylabel('Objective Function Value')
        axs[1].set_title('2nd Objective Function vs. KdRNAPAct')
        plt.tight_layout()
        plt.close()
        ret_figs.append(fig)
    
    
    if flags['auto_set_max_range']:
        rat_vals = np.linspace(min(ratios_df['actual_mRNA_ratio'].values.flatten()), max(ratios_df['actual_mRNA_ratio'].values.flatten()), 1000)

        cInh_vals = [mRNA_cActivator_to_cInhibitor(rat_val, flags['base_cActivator_val'], optimal_KdRNAPAct) for rat_val in rat_vals]
        cAct_vals = [mRNA_cInhibitor_to_cActivator(rat_val, flags['base_cInhibitor_val'], optimal_KdRNAPAct) for rat_val in rat_vals]

        try:
            flags['cActivator'] = [-4, math.log10((1+flags['additional_tolerance'])*max(cAct_vals))] # Uses a log10 range
        except:
            pass
        try:
            flags['cInhibitor'] = [-4, math.log10((1+flags['additional_tolerance'])*max(cInh_vals))] # Uses a log10 range
        except:
            pass # if there is no Inhibitor, this fails but doesn't matter
    
    return(optimal_KdRNAPAct, ret_figs)