# Modules
import numpy as np
import pandas as pd

from sympy import *
from itertools import product

# Functions
def create_grid(
        gene_exp: list,
        gene_name: list,
        equation: str,
        constants: dict,
        num_steps: int,
        **input_range: dict
    ) -> pd.DataFrame:
    """
    Create a df with the gene names, mRNA Molarity, and the corresponding grid

    :param list gene_exp: Gene expression data, NaN dropped, not log scaled
    :param list gene_name: Names corresponding to genes in gene_exp
    :param str equation: String of the form 'Eq(x,y)' that sympify can parse
    :param dict constants: Dict where the keys are symbolic variables and the 
        value is the constant's value
    :param int num_steps: Determines size of grid
    :param dict input_range: Contains names of each of the parameters as keys
        and contains log10 range of values for each parameter
    """

    # Create df with lambda functions
    lambda_df = create_lambdas(equation, constants)

    k_grid = pd.DataFrame()

    # Create a grid for every gene
    for i, expression in enumerate(gene_exp):
        # Calculate mRNA_target
        mRNA_target = np.log10(expression*10**-6*constants['mRNA_total']/constants['cell_volume']/(6.022*10**23))

        # Create list of tuples for each gene
        promoter_grid = promoter_solving(mRNA_target, num_steps, lambda_df, **input_range)

        # Create a temporary dataframe before concatenating to k_grid
        working_df = pd.DataFrame(
            columns = ['gene_name','log10_mRNA_target','parameter_grid']
        )

        working_df.at[0,'gene_name'] = gene_name[i]
        working_df.at[0,'log10_mRNA_target'] = mRNA_target
        working_df.at[0,'parameter_grid'] = promoter_grid

        k_grid = pd.concat([k_grid, working_df], ignore_index = True)
    
    k_grid = k_grid.set_index('gene_name')

    return lambda_df, k_grid

def promoter_solving(
        mRNA_target: float,
        num_steps: int,
        lambda_df: pd.DataFrame,
        **input_range: dict,
    ) -> list:
    """
    Solve for the k values

    :param dict constants: Dict where the keys are symbolic variables and the 
        value is the constant's value
    :param int num_steps: Determines size of grid
    :param pd.DataFrame lambda_df: Contains lambda functions and parameter order
    :param dict input_range: Contains names of each of the parameters as keys
        and contains log10 range of values for each parameter
    """

    # Calculate minimum k_escape and update the range
    updated_range = calculate_minimum_escape(lambda_df,mRNA_target,**input_range)

    # Create a grid from the available ranges
    working_grid = create_parameter_grid(num_steps, **updated_range)
    
    # Iterate over each combination and calculate k_eq_opening
    for i, pair in enumerate(working_grid):
        # Create dict to pass to evaluate_lambda
        values = {'mRNA': mRNA_target}
        for ii, key in enumerate(updated_range.keys()):
            values[key] = pair[ii]

        # Evaluate the lambda function
        working_grid[i] += ((evaluate_lambda('KeqOpening',lambda_df,values)),) # NOTE: UPDATE w/ Equation
    
    return np.array(object = working_grid, dtype = np.dtype([(list(input_range.keys())[0], float), (list(input_range.keys())[1], float), ('KeqOpening', float)])) # Custom dtype, might make indexing easier

def calculate_minimum_escape(
        lambda_df: pd.DataFrame,
        mRNA_target: float,
        **input_range: dict
    ) -> dict:
    """
    Calculate the minimum k_escape value and update the log_range values with it
    
    :param list kinetic_data: Generic values for kinetic equations
    :param float mRNA_target: mRNA expression data
    :param dict input_range: Contains names of each of the parameters as keys
        and contains log10 range of values for each parameter (inclusive)
    """

    # Create values dict to pass to evaluate lambda
    # NOTE: Need to update this when the equation is changed to reflect the
        # correct order in lambda_df. Run create_lambdas to see the order.
    values = {
        'KeqOpening': 2, # Assume max k_eq_opening = 100
        'KdRNAP': input_range['KdRNAP'][1], # Max k_d_RNAP
        'mRNA': mRNA_target,
        #'cActivator': input_range['cActivator'][1],
        #'KdRNAPCrp': input_range['KdRNAPCrp'][1],
    }

    min_k_escape = evaluate_lambda('kEscape',lambda_df,values)
    
    if input_range['kEscape'][1] < min_k_escape:
        raise ValueError('All values for k_escape are less than the minimum: '+
                        str(min_k_escape))
    if input_range['kEscape'][0] < min_k_escape and min_k_escape < input_range['kEscape'][1]:
        input_range['kEscape'][0] = min_k_escape
    
    return input_range

def create_parameter_grid(
        num_steps: int,
        **input_range: dict
    ) -> list:
    """
    Return a list of tuples corresponding to each combination of parameter values
    
    :param int num_steps: Determines size of grid; num_steps^(num_parameters)
    :param dict input_range: Contains names of each of the parameters as keys
        and contains log10 range of values for each parameter (inclusive)
    """

    for key in input_range:
        input_range[key] = list(np.linspace(input_range[key][0],input_range[key][1],num_steps))

    items = input_range.items()
    keys, values = zip(*items)

    return list(product(*values))

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
    values_tuple = tuple([10**values[p] for p in lambda_df.loc[solve,'order']])
    
    # Evaluate the lambda function
    return np.log10(lambda_df.loc[solve,'lambda'](values_tuple))[0]