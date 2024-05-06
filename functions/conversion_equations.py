import numpy
import pandas as pd

def get_arg_total(TFarg, TFtotal, Kdarg):
    argTot = TFarg * (-1*Kdarg + TFarg - TFtotal) / (TFarg - TFtotal)
    
    return(argTot)

def mRNA_ratio_df_to_log_tpm_df(mRNA_ratio, basal_log_tpm_val):
    return(2**(mRNA_ratio_df * basal_log_tpm_val))