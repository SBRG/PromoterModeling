# determine cActivator and cInhibior values
import pandas as pd
import numpy as np
from sympy import *
from deap import algorithms, base, creator, tools
import sys
sys.path.insert(0, '../functions/')
import promoter_solving_core as ps
import GA_core as ga
import matplotlib.pyplot as plt
import dill as pickle

def create_cAct_cInh_for_gene(ratios_df, grid_constants, eq_str, flags):
    return_figs = []
    
    # now setting RNAP concentration specific to sample
    RNAP_conc_df = pd.read_csv('../data/interim/sample_constants/RNAP_conc.csv', index_col = 0)
    f = open('../data/interim/lambda_dfs/'+flags['central_gene']+'.pkl', 'rb')
    gene_lambda_df = pickle.load(f)
    f.close()
    
    # first let's handle edge cases
    if type(flags['act_iM']) == float:
        # there is no need to run GA, solve for cInh with cAct = 0 and save those results
        #lambda_df = ps.create_lambdas(eq_str, grid_constants)
        temp = grid_constants.copy()
        temp.update({'cActivator' : flags['base_cActivator_val']})
        cInhibitors = []
        indices = []
        
        for index, row in ratios_df.iterrows():
            lambda_df = gene_lambda_df.loc[index].values[0]
            lambda_df['lambda'] = [lambdify([lambda_df.loc[index,'order']], lambda_df.loc[index,'equation']) for index in lambda_df.index]
            indices.append(index)
            if flags['use_actual_mRNA']:
                mRNA = row['actual_mRNA_ratio']
            else:
                mRNA = row['ratio_MA_inhibitor']
            temp2 = temp.copy()
            temp2.update({'mRNARatio' : mRNA})
            inputs = tuple([temp2[p] for p in lambda_df.loc['cInhibitor','order']])
            cInhibitors.append(lambda_df.loc['cInhibitor']['lambda'](inputs)[0])
        cActivators = [flags['base_cActivator_val'] for _ in cInhibitors]
        
        vals_for_GAMs = pd.DataFrame(index = indices,
                                    columns = ['cAct', 'cInh'],)

        vals_for_GAMs.cAct = list(cActivators)
        vals_for_GAMs.cInh = list(cInhibitors)

        # sanity plot
        fig = plt.figure()
        fig.suptitle('No cAct, cInh Results', fontsize = 16)
        plt.scatter(ratios_df['actual_mRNA_ratio'], vals_for_GAMs.cInh)
        plt.xlabel('actual mRNA ratio')
        plt.ylabel('cInhibitor Values')
        plt.tight_layout()
        return_figs.append(fig)
        plt.close(fig)
        
        return(return_figs, vals_for_GAMs, vals_for_GAMs)
    if type(flags['inh_iM']) == float:
        # very similar to above just flipped
        #lambda_df = ps.create_lambdas(eq_str, grid_constants)
        temp = grid_constants.copy()
        temp.update({'cInhibitor' : flags['base_cInhibitor_val']})
        cActivators = []
        indices = []
        for index, row in ratios_df.iterrows():
            lambda_df = gene_lambda_df.loc[index].values[0]
            lambda_df['lambda'] = [lambdify([lambda_df.loc[index,'order']], lambda_df.loc[index,'equation']) for index in lambda_df.index]
            indices.append(index)
            if flags['use_actual_mRNA']:
                mRNA = row['actual_mRNA_ratio']
            else:
                mRNA = row['ratio_MA_activator']
            temp2 = temp.copy()
            temp2.update({'mRNARatio' : mRNA})
            inputs = tuple([temp2[p] for p in lambda_df.loc['cActivator','order']])
            cActivators.append(lambda_df.loc['cActivator']['lambda'](inputs)[0])
        cInhibitors = [flags['base_cInhibitor_val'] for _ in cActivators]
        
        vals_for_GAMs = pd.DataFrame(index = indices,
                                    columns = ['cAct', 'cInh'],)

        vals_for_GAMs.cAct = list(cActivators)
        vals_for_GAMs.cInh = list(cInhibitors)

        # sanity plot
        fig = plt.figure()
        fig.suptitle('No cInh, cAct Results', fontsize = 16)
        plt.scatter(ratios_df['actual_mRNA_ratio'], vals_for_GAMs.cAct)
        plt.xlabel('actual mRNA ratio')
        plt.ylabel('cActivators Values')
        plt.tight_layout()
        return_figs.append(fig)
        plt.close(fig)
        
        return(return_figs, vals_for_GAMs, vals_for_GAMs)
    
    

    # setup
    rng = np.random.default_rng(seed = int(flags['seed']))
    
    # DataFrame to hold the Grid
    grid = pd.DataFrame(columns = ['mRNA_ratio','grid'], index = ratios_df.index)
    grid.loc[:,'mRNA_ratio'] = ratios_df.loc[:,'actual_mRNA_ratio']

    # Load the equation
    # NOTE: This equation was generated using Mathematica and Dan's Mathematica to Python converter function
    equation = sympify(eq_str)


    # If both exist, need to get GA going    
    cAct_range = {'cActivator': flags['cActivator']} # Use a log10 range
    cInh_range = {'cInhibitor': flags['cInhibitor']} # Use a log10 range and convert back after creating grid
    for i, condition in enumerate(grid.index):
        lambda_df = gene_lambda_df.loc[condition].values[0]
        lambda_df['lambda'] = [lambdify([lambda_df.loc[index,'order']], lambda_df.loc[index,'equation']) for index in lambda_df.index]

        
        # Create a working grid based on cActivator, we will add cInhibitor values 
        # to it to ensure they always result in mRNA ratio
        cAct_grid = ps.create_parameter_grid(num_steps = 101, **cAct_range)
        cAct_grid = [[10**x[0]] for x in cAct_grid]
        cInh_grid = ps.create_parameter_grid(num_steps = 101, **cInh_range)
        cInh_grid = [[10**x[0]] for x in cInh_grid]

        # Use a dict just in case order of tuple to sub into lambda function ever changes
        values = grid_constants.copy()
        values.update({'mRNARatio': grid.loc[condition,'mRNA_ratio']})
        for ii, pair in enumerate(cAct_grid):
            # this is untested, likely will cause issues, I'm trying to intergate the sample-specific RNAP values, but this isn't in use currently
            values['cActivator'] = pair[0] # Add cAct to values dict

            # Create a tuple in the correct order to pass into the lambda function
            values_tuple = tuple([values[p] for p in lambda_df.loc['cInhibitor','order']])

            # Evaluate the lambda function
            cAct_grid[ii] = (cAct_grid[ii][0], (lambda_df.loc['cInhibitor','lambda'](values_tuple))[0])

        for ii, pair in enumerate(cInh_grid):
            # this is untested, likely will cause issues, I'm trying to intergate the sample-specific RNAP values, but this isn't in use currently
            values['cInhibitor'] = pair[0] # Add cInh to values dict

            # Create a tuple in the correct order to pass into the lambda function
            values_tuple = tuple([values[p] for p in lambda_df.loc['cActivator','order']])

            # Evaluate the lambda function
            cInh_grid[ii] = ((lambda_df.loc['cActivator','lambda'](values_tuple))[0], cInh_grid[ii][0]) # Need to reverse the tuples to maintain (cAct, cInh) order when combining the two grids

        working_grid = sorted(cAct_grid + cInh_grid)

        # remove out of range elements
        new_work = []
        for vals in working_grid:
            cAct_val, cInh_val = vals
            if cAct_val < 10**cAct_range['cActivator'][0] or cAct_val > 10**cAct_range['cActivator'][1]:
                continue
            if cInh_val < 10**cInh_range['cInhibitor'][0] or cInh_val > 10**cInh_range['cInhibitor'][1]:
                continue
            new_work.append(vals)
        working_grid = new_work

        # Remove negative elements from working_grid
        if flags['neg_grid_toss_OR_zero'] == 'toss':
            working_grid = [(cAct, cInh) for (cAct, cInh) in working_grid if cAct >= 0 and cInh >= 0]
        elif flags['neg_grid_toss_OR_zero'] == 'zero':
            new = []
            for cAct, cInh in working_grid:
                if cAct >= 0 and cInh >= 0:
                    new.append((cAct, cInh))
                elif cAct < 0 and cInh < 0:
                    new.append((0, 0))
                elif cAct < 0:
                    new.append((0, cInh))
                else: # this is cInh is negative and cAct isn't
                    new.append((cAct, 0))
            working_grid = new

        # Save to grid df
        grid.at[condition, 'grid'] = working_grid

    # remove empty grids
    to_remove = []
    max_grid = max([len(g) for g in grid.grid])
    for index, row in grid.iterrows():
        g = row['grid']
        if len(g) == 0:
            to_remove.append(index)

    # some conditions just aren't mappable (TF KO's being one), let's drop them for now
    grid = grid.drop(to_remove)
    ratios_df = ratios_df.drop(to_remove)


    # sanity plot 1
    if flags['sanity_plots']:
        i = 0 # The integer corresponding to the condition of interest index

        x = [cAct for (cAct, _) in grid.iloc[i,1]]
        y = [cInh for (_, cInh) in grid.iloc[i,1]]

        # Start with a square Figure.
        fig = plt.figure(figsize=(3, 3))
        # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
        # the size of the marginal axes and the main axes in both directions.
        # Also adjust the subplot parameters for a square plot.
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        # Create the Axes.
        ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y, alpha=0.25)

        x_bins = np.logspace(np.log10(min(x)), np.log10(max(x)), 50)
        y_bins = np.logspace(np.log10(min(y)), np.log10(max(y)), 50)

        ax_histx.hist(x, bins=x_bins)
        ax_histy.hist(y, bins=y_bins, orientation='horizontal')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel('cAct Grid Values')
        ax.set_ylabel('cInh Grid Values')

        fig.suptitle('Condition: '+str(grid.index[i]))
        return_figs.append(fig)
        plt.close(fig)
    
    
    
    # Create fitness and individual objects
    creator.create(name = 'fitness',
                   base = base.Fitness,
                   weights = (1.0, -1.0,)) # Set to maximize Spearman correlation of MA_activator and cActivator, and minimize MA_inhibitor and cInhibitor

    creator.create(name = 'individual',
                   base = np.ndarray,
                   shape = (len(grid),), # Number of conditions
                   dtype = np.dtype([('act', float), ('inh', float)]), # Custom dtype
                   fitness = creator.fitness)

    # Import toolbox
    toolbox = base.Toolbox()

    # Register the individual and population functions
    toolbox.register(alias = 'individual',
                     function = ga.generate_individual,
                     individual_class = creator.individual,
                     grid = grid.grid,
                     rng = rng)

    toolbox.register('population',
                     tools.initRepeat,
                     list,
                     toolbox.individual)

    # Register the evaluation function
    toolbox.register(alias = 'evaluate',
                    function = ga.spearman_objective,
                    MA_df = ratios_df.loc[:,['MA_activator','MA_inhibitor']])

    # Register the selection algorithm
    toolbox.register(alias = "select", 
                     function = tools.selNSGA2, 
                     nd = 'log') 
    # I've been using selNSGA2 since it seems to run faster
    #toolbox.register("select", tools.selSPEA2)

    # Register the mutation function
    toolbox.register(alias = 'mutate', 
                     function = ga.mutate, 
                     prob = flags['mt_prob'], 
                     grid = grid.grid, 
                     rng = rng)

    # Register the crossover function
    cx_prob = 0.6 # NOTE: These values were chosen based on brute_force.ipynb
    toolbox.register(alias = "mate", 
                     function = ga.crossover, 
                     prob = flags['cx_prob'],
                     rng = rng)

    # Set the statistics to record the best individual score of each generation in 
    # the logbook
    stats = tools.Statistics(key=lambda ind: np.subtract(ind.fitness.values[0],
                                                         ind.fitness.values[1]))

    # Run the GA
    pop, logbook = ga.mu_plus_lambda(pop = toolbox.population(n = int(flags['n_ind'])), 
                                     toolbox = toolbox, 
                                     rng = rng, 
                                     mu = int(flags['mu']), 
                                     lambda_ = int(flags['lambda_']), 
                                     cxpb = flags['cxpb'], 
                                     mutpb = flags['mutpb'], 
                                     n_gen = int(flags['n_gen']), 
                                     stats = stats, 
                                     verbose = flags['verbose'])


    total_score, total_sort = ga.best_individual(pop)


    # save the non-greedy version
    GAMs_individual = pop[total_sort[-1]]

    vals_for_GAMs = pd.DataFrame(index = ratios_df.index,
                                 columns = ['cAct', 'cInh'],)

    vals_for_GAMs.cAct = list(GAMs_individual['act'])
    vals_for_GAMs.cInh = list(GAMs_individual['inh'])
    
    
    # sanity plot for GA
    if flags['sanity_plots']:
        fig, (ax1, ax2) = ga.scatter_individual(ind_one = pop[total_sort[-1]],
                                         MA = ratios_df.loc[:,['MA_activator','MA_inhibitor']],
                                         GA_parameters = flags)
        fig.suptitle('GA Results', fontsize = 16)
        plt.tight_layout()
        return_figs.append(fig)
        plt.close(fig)


    
    
    # run greedy algo
    if flags['run_greedy']:
        greedy_pop = ga.greedy_algorithm(base_individual = pop[total_sort[-1]], 
                                     n_iterations = flags['n_iter'],
                                     grid = grid.grid,
                                     toolbox = toolbox,
                                     max_steps = flags['max_steps'],
                                     n_rounds = flags['n_rounds'])

        greedy_score, greedy_sort = ga.best_individual(greedy_pop)
        greedy_score[greedy_sort[-1]]

        greedy_voting = ga.voting(population = greedy_pop,
                              grid = grid.grid)

        # Convert from condition integer index to grid tuple to create mean_ind
        mean_ind = creator.individual(greedy_pop[greedy_sort[-1]])
        for i, _ in enumerate(mean_ind):
            mean_ind[i] = grid.grid[i][int(greedy_voting._mean[i])]

        mean_ind.fitness.values = toolbox.evaluate(mean_ind)

        # save greedy
        GAMs_individual = greedy_pop[greedy_sort[-1]]

        greed_vals_for_GAMs = pd.DataFrame(index = ratios_df.index,
                                     columns = ['cAct', 'cInh'],)

        greed_vals_for_GAMs.cAct = list(GAMs_individual['act'])
        greed_vals_for_GAMs.cInh = list(GAMs_individual['inh'])
        
        
        # sanity plot for greedy
        if  flags['sanity_plots']:
            fig, (ax1, ax2) = ga.scatter_individual(ind_one = greedy_pop[greedy_sort[-1]],
                                            MA = ratios_df.loc[:,['MA_activator','MA_inhibitor']],
                                            GA_parameters = None)
            fig.suptitle('Greedy Results', fontsize = 16)
            plt.tight_layout()
            return_figs.append(fig)
            plt.close(fig)
        
        
        return(return_figs, greed_vals_for_GAMs, vals_for_GAMs)
    else:
        return(return_figs, None, vals_for_GAMs)