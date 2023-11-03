# determine cActivator and cInhibior values
import pandas as pd
import numpy as np
from sympy import *
from deap import algorithms, base, creator, tools
import sys
sys.path.insert(0, '../functions/')
import promoter_solving_core as ps
import GA_core as ga



def create_cAct_cInh_for_gene(ratios_df, grid_constants, eq_str, flags):
    # setup
    rng = np.random.default_rng(seed = flags['seed'])

    # DataFrame to hold the Grid
    grid = pd.DataFrame(columns = ['mRNA_ratio','grid'], index = ratios_df.index)
    grid.loc[:,'mRNA_ratio'] = ratios_df.loc[:,'actual_mRNA_ratio']

    # Load the equation
    # NOTE: This equation was generated using Mathematica and Dan's Mathematica to Python converter function
    equation = sympify(eq_str)


    # Create lambda functions that we can plug in to
    lambda_df = ps.create_lambdas(equation, grid_constants)

    cAct_range = {'cActivator': flags['cActivator']} # Use a log10 range
    cInh_range = {'cInhibitor': flags['cInhibitor']} # Use a log10 range and convert back after creating grid

    for i, condition in enumerate(grid.index):
        # Create a working grid based on cActivator, we will add cInhibitor values 
        # to it to ensure they always result in mRNA ratio
        cAct_grid = ps.create_parameter_grid(num_steps = 101, **cAct_range)
        cAct_grid = [[10**x[0]] for x in cAct_grid]
        cInh_grid = ps.create_parameter_grid(num_steps = 101, **cInh_range)
        cInh_grid = [[10**x[0]] for x in cInh_grid]

        # Use a dict just in case order of tuple to sub into lambda function ever changes
        values = {'mRNARatio': grid.loc[condition,'mRNA_ratio']}
        for ii, pair in enumerate(cAct_grid):
            values['cActivator'] = pair[0] # Add cAct to values dict

            # Create a tuple in the correct order to pass into the lambda function
            values_tuple = tuple([values[p] for p in lambda_df.loc['cInhibitor','order']])

            # Evaluate the lambda function
            cAct_grid[ii] = (cAct_grid[ii][0], (lambda_df.loc['cInhibitor','lambda'](values_tuple))[0])

        values = {'mRNARatio': grid.loc[condition,'mRNA_ratio']}
        for ii, pair in enumerate(cInh_grid):
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
    pop, logbook = ga.mu_plus_lambda(pop = toolbox.population(n = flags['n_ind']), 
                                     toolbox = toolbox, 
                                     rng = rng, 
                                     mu = flags['mu'], 
                                     lambda_ = flags['lambda_'], 
                                     cxpb = flags['cxpb'], 
                                     mutpb = flags['mutpb'], 
                                     n_gen = flags['n_gen'], 
                                     stats = stats, 
                                     verbose = flags['verbose'])


    total_score, total_sort = ga.best_individual(pop)


    # save the non-greedy version
    GAMs_individual = pop[total_sort[-1]]

    vals_for_GAMs = pd.DataFrame(index = ratios_df.index,
                                 columns = ['cAct', 'cInh'],)

    vals_for_GAMs.cAct = list(GAMs_individual['act'])
    vals_for_GAMs.cInh = list(GAMs_individual['inh'])
    
    
    
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
        
        return(greed_vals_for_GAMs, vals_for_GAMs)
    else:
        return(None, vals_for_GAMs)