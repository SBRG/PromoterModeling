import pandas as pd
import pickle
import numpy as np
from scipy.optimize import minimize

# these values should get overwritten in the actual code, but need to be initialized
KdRNAP = 1
KeqOpening = 1
RNAP = 1
KdRNAPCrp = 1

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


def pick_KdRNAPCrp(ratios_df, basal_constants):
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
    initial_guess = KdRNAP / 2
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
