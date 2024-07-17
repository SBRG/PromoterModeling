import numpy
import pandas as pd
import ast

def get_arg_total(TFarg, TFtotal, Kdarg):
    argTot = TFarg * (-1*Kdarg + TFarg - TFtotal) / (TFarg - TFtotal)
    
    return(argTot)

def mRNA_ratio_df_to_log_tpm_df(mRNA_ratio, flags_df, log_tpm_df):
    expression_df = pd.DataFrame(index = mRNA_ratio.index, columns = mRNA_ratio.columns)
    for column in expression_df.columns:
        basal_conds = ast.literal_eval(flags_df.loc[column]['basal_conditions'])
        basal_val = log_tpm_df.loc[column][basal_conds].values[0]
        expression_df[column] = (2**(mRNA_ratio[column] * basal_val))
    return(expression_df)

def log_tpm_df_to_mRNA_ratio_df(log_tpm_df, flags_df):
    log_tpm_df = log_tpm_df.T
    ratio_df = pd.DataFrame(index = log_tpm_df.index, columns = log_tpm_df.columns)
    for column in ratio_df.columns:
        basal_conds = ast.literal_eval(flags_df.loc[column]['basal_conditions'])
        basal_val = log_tpm_df[column].loc[basal_conds].values[0]
        ratio_df[column] = 2**log_tpm_df[column] / 2**basal_val
    return(ratio_df)