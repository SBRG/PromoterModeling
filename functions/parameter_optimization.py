import pandas as pd
import pickle
import numpy as np
from scipy.optimize import minimize
from sympy import *
from itertools import product
import math

# these values should get overwritten in the actual code, but need to be initialized
KdRNAP = None
KeqOpening = None
RNAP = None
KdRNAPCrp = None
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


# back to my code
def create_shared_lambda_df(equation_string, grid):
    global lambda_df, KdRNAP, KeqOpening, RNAP
    our_grid = grid.copy()
    equation = sympify(equation_string)
    if 'KdRNAPCrp' in our_grid:
        del(our_grid['KdRNAPCrp'])
        
    # Create lambda functions that we can plug in to
    lambda_df = create_lambdas(equation, our_grid)
    
    KdRNAP = our_grid['KdRNAP']
    KeqOpening = our_grid['KeqOpening']
    RNAP = our_grid['RNAP']

def mRNA_cInhibitor_to_cActivator(mRNA, cInhibitor, KdRNAPCrp):
    cActivator = evaluate_lambda('cActivator', lambda_df, {'cInhibitor': cInhibitor, 'mRNARatio' : mRNA, 'KdRNAPCrp' : KdRNAPCrp})
    return(cActivator)

def mRNA_cActivator_to_cInhibitor(mRNA, cActivator, KdRNAPCrp):
    cInhibitor = evaluate_lambda('cInhibitor', lambda_df, {'cActivator': cActivator, 'mRNARatio' : mRNA, 'KdRNAPCrp' : KdRNAPCrp})
    return(cInhibitor)

def cActivator_cInhibitor_to_mRNA(cActivator, cInhibitor, KdRNAPCrp):
    mRNA = evaluate_lambda('mRNARatio', lambda_df, {'cActivator': cActivator, 'cInhibitor' : cInhibitor, 'KdRNAPCrp' : KdRNAPCrp})
    return(mRNA)


def pick_KdRNAPCrp(ratios_df, flags):
    initial_guess_ratio = flags['initial_guess_ratio']
    min_cInh = flags['base_cInhibitor_val']
    # v1 assumed cInhibitor = 0, this instead sets it to a minimum value instead of zero, I'm hoping this elevates values that are below zero and is a more realistic assumption
    
    # first let's find the maximum KdRNAPCrp value for the mRNA values
    # using ratios_df is not good sampling, moving to a set spacing
    ratio_values = np.linspace(min(ratios_df['actual_mRNA_ratio']), max(ratios_df['actual_mRNA_ratio']), 1000)
    one_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - 1))
    def objective_function1(KdRNAPCrp_temp):
        
        cActs = [mRNA_cInhibitor_to_cActivator(rat_val, min_cInh, KdRNAPCrp_temp) for rat_val in ratio_values]
        
        #mRNA_to_cActivator(ratio_values, RNAP = RNAP, KdRNAP = KdRNAP, KeqOpening = KeqOpening, KdRNAPCrp = KdRNAPCrp_temp)
        

        # Set up penalty
        penalty = 0

        # I want to find the sample with the latest point above mRNA ratio 1 that has a negative cActivator value
        # i.e. I'm trying to figure out the highest KdRNAPCrp value I can use that doesn't require negatives

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
            penalty -= 1000000*KdRNAPCrp_temp
        else:
            pass#penalty += KdRNAPCrp_temp
        return(penalty)

    # setup and minimize
    initial_guess = KdRNAP * initial_guess_ratio
    bounds = [(1e-9, KdRNAP/.99)]
    result1 = minimize(objective_function1, initial_guess, bounds = bounds, tol = 1e-10)
    max_KdRNAPCrp = result1.x[0] 

    
    
    # now with our new maximum, let's look for the ideal value (and we can ignore negative value penalty I think)
    # Define your objective function
    max_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - max(ratio_values)))
    closest_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - .70*max(ratio_values)))
    def objective_function2(KdRNAPCrp_temp):
        KdRNAPCrp_temp = float(KdRNAPCrp_temp)
        cActs = [mRNA_cInhibitor_to_cActivator(rat_val, min_cInh, KdRNAPCrp_temp) for rat_val in ratio_values]

        # Set up penalty
        penalty = 0

        # let's try to maximize the 80% value
        penalty -= cActs[closest_index]

        # I want to maximize "spread" and getting a more diverse range of values... let me think on what that measure would be
        penalty -= np.corrcoef(cActs, np.arange(len(cActs)))[0, 1]
        penalty -= cActs[closest_index]

        return(penalty)

    # setup and minimize
    initial_guess = max_KdRNAPCrp / 2
    bounds = [(max_KdRNAPCrp / 100, max_KdRNAPCrp)]
    result2 = minimize(objective_function2, initial_guess, bounds = bounds, tol = 1e-10)

    # The optimal KdRNAPCrp value
    optimal_KdRNAPCrp = result2.x[0]
    
    
    if flags['auto_set_max_range']:
        rat_vals = np.linspace(min(ratios_df['actual_mRNA_ratio'].values.flatten()), max(ratios_df['actual_mRNA_ratio'].values.flatten()), 1000)

        cInh_vals = [mRNA_cActivator_to_cInhibitor(rat_val, flags['base_cActivator_val'], optimal_KdRNAPCrp) for rat_val in rat_vals]
        cAct_vals = [mRNA_cInhibitor_to_cActivator(rat_val, flags['base_cInhibitor_val'], optimal_KdRNAPCrp) for rat_val in rat_vals]

        
        flags['cActivator'] = [-2, math.log10((1+flags['additional_tolerance'])*max(cAct_vals))] # Uses a log10 range
        flags['cInhibitor'] = [-2, math.log10((1+flags['additional_tolerance'])*max(cInh_vals))] # Uses a log10 range
    
    
    return(optimal_KdRNAPCrp)




# depreceated, hopefully
if False:
    def mRNA_to_cActivator(mRNA, RNAP = RNAP, KdRNAP = KdRNAP, KdRNAPCrp = KdRNAPCrp, KeqOpening = KeqOpening): # this appears to be properly working!
        cActivator = (KdRNAPCrp*(KdRNAP + RNAP + KeqOpening*RNAP)*(-1 + \
                mRNA))/(KdRNAP*(KdRNAP + RNAP + KeqOpening*RNAP - \
                KdRNAPCrp*mRNA - RNAP*mRNA - KeqOpening*RNAP*mRNA))

        return(cActivator)

    def mRNA_to_cInhibitor(mRNA, RNAP = RNAP, KdRNAP = KdRNAP, KeqOpening = KeqOpening): # I think the mathematica code might have an error, I'm sticking with this
        cInhibitor = -(((KdRNAP + RNAP + KeqOpening*RNAP)*(-1 + mRNA))/(KdRNAP*mRNA))

        return(cInhibitor)


    def cActivator_cInhibitor_to_mRNA(cActivator, cInhibitor, RNAP = RNAP, KdRNAP = KdRNAP, KdRNAPCrp = KdRNAPCrp, KeqOpening = KeqOpening):
        mRNA = ((cActivator*KdRNAP + KdRNAPCrp)*(KdRNAP + RNAP + \
                KeqOpening*RNAP))/((1 + cActivator + cInhibitor)*KdRNAP*KdRNAPCrp + \
                cActivator*KdRNAP*(1 + KeqOpening)*RNAP + KdRNAPCrp*(1 + \
                KeqOpening)*RNAP)

        return(mRNA)

    def log_tpm_to_mRNA_conc(val):
        mRNA = (2**val)*(10**-6)*1800/cell_volume/(6.022*(10**23))
        return(mRNA)
    
    
    def pick_KdRNAPCrp(ratios_df, basal_constants, initial_guess_ratio = 0.25):
        KdRNAP = basal_constants['KdRNAP']
        KeqOpening = basal_constants['KeqOpening']
        RNAP = basal_constants['RNAP']

        # first let's find the maximum KdRNAPCrp value for the mRNA values
        # using ratios_df is not good sampling, moving to a set spacing
        ratio_values = np.linspace(min(ratios_df['actual_mRNA_ratio']), max(ratios_df['actual_mRNA_ratio']), 1000)
        one_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - 1))
        def objective_function1(KdRNAPCrp_temp):
            cActs = mRNA_to_cActivator(ratio_values, RNAP = RNAP, KdRNAP = KdRNAP, KeqOpening = KeqOpening, KdRNAPCrp = KdRNAPCrp_temp)


            # Set up penalty
            penalty = 0

            # I want to find the sample with the latest point above mRNA ratio 1 that has a negative cActivator value
            # i.e. I'm trying to figure out the highest KdRNAPCrp value I can use that doesn't require negatives

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
                penalty -= 1000000*KdRNAPCrp_temp
            else:
                pass#penalty += KdRNAPCrp_temp
            return(penalty)

        # setup and minimize
        initial_guess = KdRNAP * initial_guess_ratio
        bounds = [(1e-9, KdRNAP/.99)]
        result1 = minimize(objective_function1, initial_guess, bounds = bounds, tol = 1e-10)
        max_KdRNAPCrp = result1.x[0] 



        # now with our new maximum, let's look for the ideal value (and we can ignore negative value penalty I think)
        # Define your objective function
        max_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - max(ratio_values)))
        closest_index = min(range(len(ratio_values)), key=lambda i: abs(ratio_values[i] - .70*max(ratio_values)))
        def objective_function2(KdRNAPCrp_temp):
            cActs = mRNA_to_cActivator(ratio_values, RNAP = RNAP, KdRNAP = KdRNAP, KeqOpening = KeqOpening, KdRNAPCrp = KdRNAPCrp_temp)

            # Set up penalty
            penalty = 0

            # let's try to maximize the 80% value
            penalty -= cActs[closest_index]

            # I want to maximize "spread" and getting a more diverse range of values... let me think on what that measure would be
            penalty -= np.corrcoef(cActs, np.arange(len(cActs)))[0, 1]
            penalty -= cActs[closest_index]

            return(penalty)

        # setup and minimize
        initial_guess = max_KdRNAPCrp / 2
        bounds = [(max_KdRNAPCrp / 100, max_KdRNAPCrp)]
        result2 = minimize(objective_function2, initial_guess, bounds = bounds, tol = 1e-10)

        # The optimal KdRNAPCrp value
        optimal_KdRNAPCrp = result2.x[0]
        return(optimal_KdRNAPCrp)
