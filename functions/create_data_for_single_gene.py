"""
Contains a function which creates all necessary data files for GAMS for a specific gene

Functions:
create_data_for_gene - Inputs flags to generate all necessary data files to run GAMS for a specific gene

"""

# imports
import os
import dill as pickle
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from sympy import *
import sys
sys.path.insert(0, '../functions/')
import basal_model_calcs as bmc
import mRNA_ratios as mr
import parameter_optimization as po
import create_cAct_cInh_vals as cv
import interface_GAMS as iG


def create_data_for_gene(flags):
    """
    Inputs flags to generate all necessary data files to run GAMS for a specific gene. Generates mRNA ratios, lambda_dfs, calculates KdRNAPActivator, creates greedy and standard cActivator and cInhibitor values, saves off to various data folders
    
    Inputs:
        flags (dict) : dictionary of settings flags and constants values
    
    Returns:
        gene_figs (array) : set of figures used as sanity check on production of necessary values for GAMS
    """
        
    # setup
    gene_figs = []
    eq_str = flags['eq_str']
    if 'hard_set_save_folder' in flags:
        save_folder = flags['hard_set_save_folder']
    else:
        save_folder = '../data/output/saved_gene_results/'+flags['central_gene']
    
    ############################################################
    # create mRNA ratios and MA values
    ############################################################
    df_name = flags['central_gene']+'_zerod'+str(flags['use_zerod_A_matrix'])+'_mRNA_ratios_and_MA_vals.csv'
    if not flags['force_rerun'] and os.path.exists('../data/saved_mRNA_ratios_MA_vals/'+df_name):
        ratios_df = pd.read_csv('../data/processed/saved_mRNA_ratios_MA_vals/'+df_name, index_col = 0)
    else:
        ratios_df = mr.calculate_mRNA_ratios_and_MA_values(flags['act_iM'], flags['inh_iM'], flags)
        ratios_df.to_csv('../data/processed/saved_mRNA_ratios_MA_vals/'+df_name)
    if flags['sanity_plots']:
        # sanity check plot
        fig, axs = plt.subplots(1, 3, figsize = (10, 3))
        axs[0].hist(ratios_df.actual_mRNA_ratio)
        axs[0].set_title('mRNA ratio')
        axs[0].axvline(x = 1, c = 'k', ls = '--')
        axs[1].hist(ratios_df.MA_activator)
        axs[1].set_title('MA activator')
        axs[1].axvline(x = 0, c = 'k', ls = '--')
        axs[2].hist(ratios_df.MA_inhibitor)
        axs[2].set_title('MA_inhibitor')
        axs[2].axvline(x = 0, c = 'k', ls = '--')
        # add a big axes, hide frame
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.ylabel('Count')
        gene_figs.append(fig)
        plt.close(fig)
    if flags['only_create_ratios']:
        # save off results
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)
        os.mkdir(save_folder)

        # let's save off the relevant things in pickle
        pickle_out = open(save_folder+'/figures.pkl', 'wb')
        pickle.dump(gene_figs, pickle_out)
        pickle_out.close()

        pickle_out = open(save_folder+'/ratios_df.pkl', 'wb')
        pickle.dump(ratios_df, pickle_out)
        pickle_out.close()
        
        return(gene_figs)
    
    ############################################################
    # create lambda dfs for various calculations --- very slow
    ############################################################
    lambda_df_path = '../data/interim/lambda_dfs/'+flags['central_gene']+'.pkl'
    if flags['force_basal_rerun'] or not os.path.exists(lambda_df_path):
        # load in
        if os.path.exists(lambda_df_path):
            lambdas_df = pd.read_pickle(lambda_df_path)
        else:
            lambdas_df = pd.DataFrame(index = ratios_df.index)
        RNAP_conc_df = pd.read_csv('../data/interim/sample_constants/RNAP_conc.csv', index_col = 0)
        
        # calculate gene specific grid constants
        num_grid_steps = 3 
        if 'grid_steps' in flags:
            num_grid_steps = flags['grid_steps']
        gene_grid_constants = bmc.basal_values(eq_str, flags, num_steps = num_grid_steps)
        new_col = []
        for sample in ratios_df.index:
            # adjust constants
            gene_grid_constants.update({'RNAP' : RNAP_conc_df.loc[sample]['RNAP']})
            
            # create lambdas_df and put it in dataframe
            lambda_df = po.create_lambdas(eq_str, gene_grid_constants)
            
            # add lambda func
            for parameter in lambda_df.index:
                lambda_df.loc[parameter,'lambda'] = lambdify([lambda_df.loc[parameter,'order']],
                                                    lambda_df.loc[parameter,'equation'])
            
            new_col.append(lambda_df[['equation', 'order', 'lambda']])
        lambdas_df[flags['central_gene']] = new_col
        
        # save it
        f = open(lambda_df_path, 'wb')
        pickle.dump(lambdas_df, f)
        f.close()
    
    ############################################################
    # pick KdRNAPActivator value, limit cActivator and cInhibitor based on it
    ############################################################
    # load in calculator
    gene_grid_name = '../data/interim/gene_grid_constants/'+flags['central_gene']+'.pkl'
    if flags['only_check_KdRNAPAct'] or flags['force_rerun'] or not os.path.exists(gene_grid_name):  
        # basal model calculations
        num_grid_steps = 3
        if 'grid_steps' in flags:
            num_grid_steps = flags['grid_steps']
        grid_constants = bmc.basal_values(eq_str, flags, num_steps = num_grid_steps)

        # pick KdRNAPAct
        po.create_shared_lambda_df(eq_str, grid_constants)
        grid_constants['KdRNAPAct'], ret_figs = po.pick_KdRNAPActivator(ratios_df, flags)
        for temp in ret_figs:
            gene_figs.append(temp)

        # save off grid constants
        pickle_out = open(gene_grid_name, 'wb')
        pickle.dump(grid_constants, pickle_out)
        pickle_out.close()
    else:
        pickle_in = open(gene_grid_name, 'rb')
        grid_constants = pickle.load(pickle_in)
        pickle_in.close()
    if flags['sanity_plots']:
        # sanity check plot
        
        # loading / setup
        po.create_shared_lambda_df(eq_str, grid_constants)

        # if you get weird results here, look at egulonML/parameter_optimization/0_framework.ipynb
        # it does the same thing as the function with plots along the way
        # try adjusting the initial guess for the first optimization
        # that is line 74 of functions/parameter_optimization.py

        # however, it is a sanity check to see if these values are near-correct
        rat_vals = np.linspace(min(ratios_df['actual_mRNA_ratio'].values.flatten()), max(ratios_df['actual_mRNA_ratio'].values.flatten()), 1000)
        # calculate cInh and cAct
        cInh_vals = [po.mRNA_cActivator_to_cInhibitor(rat_val, flags['base_cActivator_val'], grid_constants['KdRNAPAct']) for rat_val in rat_vals]
        cAct_vals = [po.mRNA_cInhibitor_to_cActivator(rat_val, flags['base_cInhibitor_val'], grid_constants['KdRNAPAct']) for rat_val in rat_vals]

        fig, axs = plt.subplots(1, 2, figsize = (8, 3))
        ax1 = axs[0]
        l1, = ax1.plot(rat_vals, cInh_vals)
        plt.xlabel('mRNA Ratio')
        ax1.set_ylabel('cInhibitor', color = 'blue')
        ax1.tick_params(axis = 'y', labelcolor = 'blue')
        ax2 = ax1.twinx()
        l2, = ax2.plot(rat_vals, cAct_vals, color = 'red')
        ax2.set_ylabel('cActivator', color = 'red')
        ax2.tick_params(axis = 'y', labelcolor = 'red')
        ax1.axhline(y = 0, ls = '--', c = 'k')
        ax1.axvline(x = 1, ls = '--', c = 'k')
        # let's rescale cInhibitor (ax1) so that 0 is at the same point
        m1, M1 = ax1.get_ylim()
        percent1_up = (0 - m1) / (M1 - m1)
        m2, M2 = ax2.get_ylim()
        percent2_up = (0 - m2) / (M2 - m2)

        if percent1_up < percent2_up:
            # zero is higher than it should, so adjust it down by lowering the min
            m1 = percent2_up * M1 / (percent2_up - 1)
            ax1.set_ylim(m1, M1)
        else:
            # zero is lower than it should, so adjust it up
            M1 = m1 - (m1 / percent2_up)
            ax1.set_ylim(m1, M1)

        ax1.set_title('cAct and cInh ranges assuming other = 0')

        if flags['auto_set_max_range']:
            try:
                flags['cActivator'] = [-2, math.log10((1+flags['additional_tolerance'])*max(cAct_vals))] # Uses a log10 range
            except:
                pass
            try:
                flags['cInhibitor'] = [-2, math.log10((1+flags['additional_tolerance'])*max(cInh_vals))] # Uses a log10 range
            except:
                pass
            
        # let's create a 2D heatmap version of this, colored by the mRNA ratio
        cInh_range = np.linspace(0, max(cInh_vals), 100)
        cAct_range = np.linspace(0, max(cAct_vals), 100)
        mRNA_vals = pd.DataFrame(index = cInh_range, columns = cAct_range)
        for cInh in mRNA_vals.index:
            for cAct in mRNA_vals.columns:
                mRNA_vals.loc[cInh][cAct] = po.cActivator_cInhibitor_to_mRNA(cAct, cInh, grid_constants['KdRNAPAct'])
        mRNA_vals = mRNA_vals.T.astype(float)

        # Convert the cInh_range and cAct_range to meshgrids for plotting
        cInh, cAct = np.meshgrid(cInh_range, cAct_range)

        # Create the heatmap
        try:
            heatmap = axs[1].pcolormesh(cInh, cAct, mRNA_vals.values, shading='auto', cmap='viridis')
            plt.colorbar(heatmap, label='mRNA values')
            axs[1].set_xlabel('cInhibitor')
            axs[1].set_ylabel('cActivator')
            axs[1].set_title('2D Heatmap of mRNA values')
            plt.tight_layout()
            gene_figs.append(fig)
            plt.close(fig)
        except:
            pass
            # this means something went wrong in creating the range, presumably because cInh or cAct contains non-finite values
    if flags['only_check_KdRNAPAct']:
        # save off results
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)
        os.mkdir(save_folder)

        # save off constants used
        pickle_out = open(save_folder+'/constants.pkl', 'wb')
        pickle.dump(grid_constants, pickle_out)
        pickle_out.close()

        # let's save off the relevant things in pickle
        pickle_out = open(save_folder+'/figures.pkl', 'wb')
        pickle.dump(gene_figs, pickle_out)
        pickle_out.close()

        pickle_out = open(save_folder+'/ratios_df.pkl', 'wb')
        pickle.dump(ratios_df, pickle_out)
        pickle_out.close()
        
        return(gene_figs)
    
    ############################################################
    # determine cActivator and cInhibior values, and greedy
    ############################################################
    greedy_path = '../data/processed/cAct_cInh_vals/'+flags['central_gene']+'_greedy.pkl'
    norm_path = '../data/processed/cAct_cInh_vals/'+flags['central_gene']+'.pkl'
    return_figs = []
    if flags['force_rerun']:
        return_figs, greedy_cAct_cInh_df, cAct_cInh_df = cv.create_cAct_cInh_for_gene(ratios_df, grid_constants, flags)
        if flags['run_greedy']:
            pickle_out = open(greedy_path, 'wb')
            pickle.dump(greedy_cAct_cInh_df, pickle_out)
            pickle_out.close()
            pickle_out = open(norm_path, 'wb')
            pickle.dump(cAct_cInh_df, pickle_out)
            pickle_out.close()
        else:
            pickle_out = open(norm_path, 'wb')
            pickle.dump(cAct_cInh_df, pickle_out)
            pickle_out.close()
    else:
        if flags['run_greedy'] and os.path.exists(norm_path):
            pickle_in = open(norm_path, 'rb')
            cAct_cInh_df = pickle.load(pickle_in)
            pickle_in.close()
        elif os.path.exists(norm_path) and os.path.exists(greedy_path):
            pickle_in = open(norm_path, 'rb')
            cAct_cInh_df = pickle.load(pickle_in)
            pickle_in.close()
            pickle_in = open(greedy_path, 'rb')
            greedy_cAct_cInh_df = pickle.load(pickle_in)
            pickle_in.close()
        else: # need to rerun
            return_figs, greedy_cAct_cInh_df, cAct_cInh_df = cv.create_cAct_cInh_for_gene(ratios_df, grid_constants, flags)
            if flags['run_greedy']:
                pickle_out = open(greedy_path, 'wb')
                pickle.dump(greedy_cAct_cInh_df, pickle_out)
                pickle_out.close()
                pickle_out = open(norm_path, 'wb')
                pickle.dump(cAct_cInh_df, pickle_out)
                pickle_out.close()
            else:
                pickle_out = open(norm_path, 'wb')
                pickle.dump(cAct_cInh_df, pickle_out)
                pickle_out.close()
    for fig in return_figs:
        gene_figs.append(fig)
    
    ############################################################
    # save off results
    ############################################################
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    os.mkdir(save_folder)

    # save off constants used
    pickle_out = open(save_folder+'/constants.pkl', 'wb')
    pickle.dump(grid_constants, pickle_out)
    pickle_out.close()

    # let's save off the relevant things in pickle
    pickle_out = open(save_folder+'/figures.pkl', 'wb')
    pickle.dump(gene_figs, pickle_out)
    pickle_out.close()

    pickle_out = open(save_folder+'/ratios_df.pkl', 'wb')
    pickle.dump(ratios_df, pickle_out)
    pickle_out.close()

    pickle_out = open(save_folder+'/cAct_cInh.pkl', 'wb')
    pickle.dump(cAct_cInh_df, pickle_out)
    pickle_out.close()

    if flags['run_greedy']:
        try:
            pickle_out = open(save_folder+'/greedy_cAct_cInh.pkl', 'wb')
            pickle.dump(greedy_cAct_cInh_df, pickle_out)
            pickle_out.close()
        except:
            pass # Greedy sometimes not run, just don't save if the case
    
    
    return(gene_figs)