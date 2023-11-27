import numpy
import pandas as pd

def get_arg_total(TFarg, TFtotal, Kdarg):
    argTot = TFarg * (-1*Kdarg + TFarg - TFtotal) / (TFarg - TFtotal)
    
    return(argTot)