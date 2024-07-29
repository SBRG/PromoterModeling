# Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from deap import base, tools

# NOTE: Below are the functions used by the deap toolbox. The creator and 
# toolbox must be registered in each notebook they are used in.

def generate_individual(
        individual_class: type, 
        grid: pd.Series, 
        rng: object
        ) -> object:
    '''
    Randomly pick tuples from each condition's grid to create an individual
    
    :param type individual_class: Class that the individual will inherit from
    :param pd.Series grid: Contains the parameter grids for each condition
    :param object rng: np.random.generator or np.random.RandomState class
    '''
    individual = np.empty(shape=individual_class.shape, 
                          dtype=individual_class.dtype)

    for i, condition in enumerate(grid):
        individual[i] = tuple(rng.choice(a=condition, 
                                         size=1, 
                                         replace=False)[0])
        
    return individual_class(individual)

def spearman_objective(
        individual: object, 
        MA_df: pd.DataFrame
        ) -> tuple:
    """
    Calculate spearman coefficient between MA_act and cAct as well as MA_inh and
        cInh
    
    :param object individual: DEAP individual
    :param pd.DataFrame MA_df: Df with columns = ['MA_activator','MA_inhibitor']
    """

    # Create arrays with the ordered condition parameters
    MA_activator = MA_df.loc[:,'MA_activator']
    MA_inhibitor = MA_df.loc[:,'MA_inhibitor']

    ind_activator = individual['act']
    ind_inhibitor = individual['inh']
    
    # Calculate the spearman rank coefficient
    activator_spearman = spearmanr(MA_activator, ind_activator)[0]
    inhibitor_spearman = spearmanr(MA_inhibitor, ind_inhibitor)[0]
    
    return activator_spearman, inhibitor_spearman,

def mutate(
        individual: object, 
        prob: float, 
        grid: pd.Series, 
        rng: object
        ) -> object:
    '''
    Given that the individual is chosen to undergo mutation, go through each
        condition in the individual and replace the tuple with a randomly chosen
        tuple if the rng result is less than the set probability
    
    :param object individual: Class that the individual will inherit from
    :param float prob: Probability that a given condition will be mutated
    :param pd.Series grid: Contains the parameter grids for each condition
    :param object rng: np.random.generator or np.random.RandomState class
    '''

    # Iterate over all conditions if an individual is selected to mutate
    for i, _ in enumerate(individual):
        if rng.random() < prob:
            # Select a new set of parameters for a condition from the grid
            # Requires the "grid" to be in a specific location of the dataframe
            individual[i] = tuple(rng.choice(a=grid.iloc[i], 
                                             size=1, 
                                             replace=False)[0])

    return individual

def crossover(
        ind_one: object, 
        ind_two: object, 
        prob: float, 
        rng: object
        ) -> list:
    '''
    Given two randomly selected individuals, swap the tuples of a condition if 
        the rng result is less than the set probability
    
    :param object ind_one: The first individual undergoing crossover
    :param object ind_two: The second individual undergoing crossover
    :param float prob: Probability that a given condition will be mutated
    :param object rng: np.random.generator or np.random.RandomState class
    '''

    # Individuals are already deep-copied before going into this function
    for i, _ in enumerate(ind_one):
        if rng.random() < prob:
            # Use copies to avoid modifying the np.array in place
            ind_one[i], ind_two[i] = ind_two[i].copy(), ind_one[i].copy()

    return ind_one, ind_two

def best_individual(
        pop: list
        ) -> list:
    """
    Returns a list containing the total scores for each individual and a list 
    containining the np.argsort indices for the ascending sorted scores

    :param list pop: List of individuals with fitness scores
    """
    a, b = zip(*[pop[i].fitness.values for i in range(len(pop))])
    total_scores = np.subtract(a,b) # Since objective weights are (1.0, -1.0)
    sorted_index = np.argsort(total_scores)

    return total_scores, sorted_index

def var_or(
        pop: list, 
        toolbox: object, 
        lambda_: int, 
        cxpb: float, 
        mutpb: float,
        rng: object
        ) -> list:
    '''
    NOTE: This function was copied from DEAP and modified to work with the
        np.random.generator

    Part of an evolutionary algorithm applying only the variation part
    (crossover, mutation **or** reproduction). The modified individuals have
    their fitness invalidated. The individuals are cloned so returned
    pop is independent of the input pop.

    :param list pop: A list of individuals to vary.
    :param object toolbox: A :class:`~deap.base.Toolbox` that contains the 
        evolution operators.
    :param int lambda\_: The number of children to produce
    :param float cxpb: The probability of mating two individuals.
    :param float mutpb: The probability of mutating an individual.
    :param object rng: np.random.generator or np.random.RandomState class

    The variation goes as follow. On each of the *lambda_* iteration, it
    selects one of the three operations; crossover, mutation or reproduction.
    In the case of a crossover, two individuals are selected at random from
    the parental pop :math:`P_\mathrm{p}`, those individuals are cloned
    using the :meth:`toolbox.clone` method and then mated using the
    :meth:`toolbox.mate` method. Only the first child is appended to the
    offspring pop :math:`P_\mathrm{o}`, the second child is discarded.
    In the case of a mutation, one individual is selected at random from
    :math:`P_\mathrm{p}`, it is cloned and then mutated using using the
    :meth:`toolbox.mutate` method. The resulting mutant is appended to
    :math:`P_\mathrm{o}`. In the case of a reproduction, one individual is
    selected at random from :math:`P_\mathrm{p}`, cloned and appended to
    :math:`P_\mathrm{o}`.

    This variation is named *Or* because an offspring will never result from
    both operations crossover and mutation. The sum of both probabilities
    shall be in :math:`[0, 1]`, the reproduction probability is
    1 - *cxpb* - *mutpb*.
    '''

    assert (cxpb + mutpb) <= 1.0, (
        "The sum of the crossover and mutation probabilities must be smaller "
        "or equal to 1.0.")

    offspring = []
    for _ in range(lambda_):
        op_choice = rng.random()
        if op_choice < cxpb:
            # Apply crossover
            # Use rng.choice with arange to get a random integer and copy that 
            # individual, rather than generate an np.array and create an 
            # individual from that
            ind1, ind2 = [toolbox.clone(pop[i]) for i in rng.choice(a=np.arange(start=0, 
                                                                                stop=len(pop), 
                                                                                step=1), 
                                                                    size=2, 
                                                                    replace=False)]
            ind1, ind2 = toolbox.mate(ind1, ind2)
            del ind1.fitness.values
            del ind2.fitness.values
            offspring.append(ind1) 
            # We currently only take one of the crossover indidivuals
        elif op_choice < cxpb + mutpb:
            # Apply mutation
            ind = toolbox.clone(pop[rng.choice(a=np.arange(start=0, 
                                                           stop=len(pop), 
                                                           step=1), 
                                               size=1, 
                                               replace=False)[0]])
            ind = toolbox.mutate(ind)
            del ind.fitness.values
            offspring.append(ind)
        else:
            # Apply reproduction
            offspring.append(rng.choice(a=pop, 
                                        size=1, 
                                        replace=False))

    return offspring

def mu_plus_lambda(pop = list,
                   toolbox = base.Toolbox(),
                   rng = object,
                   mu = int,
                   lambda_ = int,
                   cxpb = float,
                   mutpb = float,
                   n_gen = int,
                   n_iter = None,
                   grid = None,
                   stats = None,
                   hall_of_fame = None,
                   verbose = __debug__):
    """
    Modified DEAP mu+lambda evolutionary algorithm using var_or

    :param list pop: List of individuals to serve as the starting pop
    :param base.Toolbox() toolbox: DEAP class containing evolution operators
    :param object rng: np.random.generator or np.random.RandomState class
    :param int mu: Number of individuals to select for the next generation
    :param int lambda_: Number of children to produce at each generation
    :param float cxpb: Probability that an offspring is produced by crossover
    :param float mutpb: Probability that an offspring is produced by mutation
    :param int n_gen: Number of generations to run
    :param int n_iter: Greedy offspring pop size
    :param pd.Series grid: Series containing grids for each condition
    :param stats: DEAP class containing the types of statistics to record in the 
        logbook
    :param halloffame: DEAP class containing the best individuals evaluated
    :param verbose: Whether or not to print statistics for each generation
    :returns list pop: Final pop
    :returns logbook: DEAP class containing stats for every generation

    evaluate(pop)
    for g in range(ngen):
        offspring = varOr(pop, toolbox, lamda_, cxpb, mutpb)
        evaluate(offspring)
        gradient_offspring = create_gradient_offspring(hof[0], toolbox)
        evaluate(gradient_offspring)
        pop = select(pop+offspring+gradient_offspring, mu)
    """

    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals', 'best'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in pop if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, 
                            invalid_ind)
    for ind, fit in zip(invalid_ind, 
                        fitnesses):
        ind.fitness.values = fit

    if hall_of_fame is not None:
        hall_of_fame.update(pop)

    total_scores, sorted_index = best_individual(pop)

    record = stats.compile(pop) if stats is not None else {}
    logbook.record(gen=0, 
                   nevals=len(invalid_ind), 
                   best=total_scores[sorted_index[-1]], 
                   **record)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(1, n_gen + 1):
        # Vary the pop
        offspring = var_or(pop, 
                           toolbox, 
                           lambda_, 
                           cxpb, 
                           mutpb,
                           rng)

        '''
        # NOTE: No longer applying in this fashion, apply after GA has finished
        # Greedy offspring
        if gen % 100 == 0:
            greedy_offspring, _ = create_greedy_offspring(pop[sorted_index[-1]], 
                                                          n_iter, 
                                                          grid)
            offspring = offspring+greedy_offspring
        '''

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, 
                                invalid_ind)
        for ind, fit in zip(invalid_ind, 
                            fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        # NOTE: I feel like HoF should update after the next pop is created? 
        # That way we only have to compare max mu individuals to the HoF?
        if hall_of_fame is not None:
            hall_of_fame.update(offspring)

        # Manually ensure elitism
        pop = pop + offspring
        _, temp_sorted = best_individual(pop)
        hof_individual = pop[temp_sorted[-1]]

        # Select the next generation pop
        pop[:] = toolbox.select(pop, 
                                mu)
        pop.append(hof_individual)

        # Update the statistics with the new pop
        total_scores, sorted_index = best_individual(pop)

        record = stats.compile(pop) if stats is not None else {}
        logbook.record(gen=gen, 
                       nevals=len(invalid_ind), 
                       best=total_scores[sorted_index[-1]], 
                       **record)
        if verbose:
            print(logbook.stream)

    return pop, logbook

def greedy_algorithm(base_individual: object,
                     n_iterations: int,
                     grid: pd.Series,
                     toolbox = base.Toolbox(),
                     max_steps: int = 30,
                     n_rounds: int = 100):
    """
    Returns a population of modified individuals that have different parameters
        for one condition

    :param object base_individual: Individual to copy parameters from
    :param int n_iterations: How many random copies should be created
    :param pd.Series grid: Series containing grids for each condition
    :param base.Toolbox() toolbox: DEAP class containing evolution operators
    :param int max_steps: Maximum number of steps the greedy algorithm can take
        at a time before moving on to a different condition
    """
    n_iterations = int(n_iterations)
    max_steps = int(max_steps)
    n_rounds = int(n_rounds)
    
    # Create population to hold the individuals we are trying out
    population = []

    # Create a pd.Series to hold position of each condition's parameter tuple
    position = pd.Series(index = grid.index, dtype = 'Int64')

    for i, (ind_act, ind_inh) in enumerate(base_individual):
        # Iterate through the grid parameters for each condition
        for ii, (grid_act, grid_inh) in enumerate(grid[i]):
            if (ind_act, ind_inh) == (grid_act, grid_inh):
                position[i] = ii

    # Record the number of steps needed for each condition at each iteration
    steps = pd.DataFrame(index=grid.index, columns=range(n_iterations))

    # Shuffle the order of the conditions and run it many times
    for i in range(n_iterations):
        # Create the iteration's temporary individual
        base_score = np.subtract(base_individual.fitness.values[0], 
                                 base_individual.fitness.values[1])
        temp_individual = toolbox.clone(base_individual)

        # Baseline score
        old_l_score = base_score
        old_r_score = base_score

        # Create a copy of the grid and position series as new objects
        temp_position = position.copy(deep=True)

        for ii in range(n_rounds):
            # Copy the grid and position series and shuffle them
            temp_grid = grid.sample(frac = 1, 
                                    replace = False, 
                                    random_state = (i+1)*(ii+1)) # NOTE: What's the best way to set the random state?

            # Loop over the shuffled conditions
            for cond in temp_grid.index:
                # Create the condition's temporary individuals
                l_individual = toolbox.clone(temp_individual)
                r_individual = toolbox.clone(temp_individual)

                # Identify where in the individual the condition's parameters are
                cond_loc = int(grid.index.get_loc(cond))

                if int(temp_position.loc[cond]) == 0:
                    l_position = 0
                else:
                    l_position = int(temp_position.loc[cond]-1)
                    l_steps = 1

                    # Search the "left" side
                    while l_position > -1:
                    # Modify the individual, evaluate, and score
                        l_individual[cond_loc] = temp_grid.loc[cond][l_position]
                        l_individual.fitness.values = toolbox.evaluate(l_individual)
                        new_l_score = np.subtract(l_individual.fitness.values[0],
                                                  l_individual.fitness.values[1])

                        # Compare the score against the best individual
                        if new_l_score > old_l_score:
                            # Continue if not at the end of the grid
                            if l_position == 0 or l_steps == max_steps:
                                old_l_score = new_l_score
                                break
                            else:
                                l_position -= 1
                                l_steps += 1
                        else: # Roll the individual back one step
                            l_position += 1
                            #print(i, ii, cond, l_position, len(temp_grid.loc[cond]))
                            l_individual[cond_loc] = temp_grid.loc[cond][l_position]
                            l_individual.fitness.values = toolbox.evaluate(l_individual)
                            old_l_score = np.subtract(l_individual.fitness.values[0],
                                                      l_individual.fitness.values[1])
                            break

                if int(temp_position.loc[cond])+1 == len(temp_grid.loc[cond]):
                    # This if statement is needed to prevent it from adding one, realizing it is out of bounds, performing better than the left side, and updating to a grid location that is out of bounds
                    r_position = int(temp_position.loc[cond])
                else:
                    r_position = int(temp_position.loc[cond]+1)
                    r_steps = 1
                    
                    # Search the "right" side
                    while r_position < len(temp_grid.loc[cond]):
                        # Modify the individual, evaluate, and score
                        r_individual[cond_loc] = temp_grid.loc[cond][r_position]
                        r_individual.fitness.values = toolbox.evaluate(r_individual)
                        new_r_score = np.subtract(r_individual.fitness.values[0],
                                                  r_individual.fitness.values[1])

                        # Compare the score against the best individual
                        if new_r_score > old_r_score:
                            # Continue if not at the end of the grid
                            if r_position+1 == len(temp_grid.loc[cond]) or r_steps == max_steps:
                                old_r_score = new_r_score
                                break
                            else:
                                r_position += 1
                                r_steps += 1
                        else: # Roll the individual back one step
                            r_position -= 1
                            r_individual[cond_loc] = temp_grid.loc[cond][r_position]
                            r_individual.fitness.values = toolbox.evaluate(r_individual)
                            old_r_score = np.subtract(r_individual.fitness.values[0],
                                                      r_individual.fitness.values[1])
                            break

                # Compare r_individual and l_individual, use best, update the positions for the next round
                if old_r_score > old_l_score:
                    temp_individual = r_individual
                    temp_position.loc[cond] = r_position
                else:
                    temp_individual = l_individual
                    temp_position.loc[cond] = l_position
                temp_individual.fitness.values = toolbox.evaluate(temp_individual)
            
        # Add individual to population
        population.append(temp_individual)
    
    return population

def voting(population: list, 
           grid: pd.Series,
           ) -> pd.DataFrame:
    '''
    Create a pd.DataFrame containing the rounded mean, variance, and rounded 
        median of each condition's integer index corresponding to the grid value 
    '''

    # Create a pd.Series to hold position of each condition's parameter tuple
    position = pd.DataFrame(index=grid.index, columns=range(len(population)), dtype='Int64')

    for i, ind in enumerate(population):
        for ii, (ind_act, ind_inh) in enumerate(ind):
            # Iterate through the grid parameters for each condition
            for iii, (grid_act, grid_inh) in enumerate(grid[ii]):
                if (ind_act, ind_inh) == (grid_act, grid_inh):
                    position.iloc[ii, i] = iii

    position['_mean'] = position.iloc[:,0:len(population)].mean(axis=1).round()
    position['_var'] = position.iloc[:,0:len(population)].var(axis=1)
    position['_med'] = position.iloc[:,0:len(population)].median(axis=1).round()
    
    return position

# Functions related to plotting results
def scatter_individual(ind_one: object, 
                       ind_two: object = None, 
                       MA = pd.DataFrame,
                       GA_parameters: dict = None):
    '''
    Create scatter plots of cAct vs MA_activator and cInh vs MA_inhibitor

    :param object ind_one: The individual to plot
    :param None ind_two: Second individual to plot on the same axes
    :param pd.DataFrame MA: DataFrame containing MA_activator and MA_inhibitor
    :param dict GA_parameters: Dictionary containing parameters to include in 
        the name of the plot
    '''

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), squeeze=True)

    if GA_parameters is not None:
        title = 'seed='+str(GA_parameters['seed'])+', n_gen='+str(GA_parameters['n_gen'])+', pop='+str(GA_parameters['n_ind'])+', μ='+str(GA_parameters['mu'])+', λ='+str(GA_parameters['lambda_'])+', cxpb='+str(GA_parameters['cxpb'])+', mutpb='+str(GA_parameters['mutpb'])+', mt_prob='+str(GA_parameters['mt_prob'])+', cx_prob='+str(GA_parameters['cx_prob'])
        fig.suptitle(title)

    ax1.set_xlabel('MA_activator')
    ax2.set_xlabel('MA_inhibitor')
    ax1.set_ylabel('cActivator')
    ax2.set_ylabel('cInhibitor')
    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax1.set_title('Score: '+str(ind_one.fitness.values[0]))
    ax2.set_title('Score: '+str(ind_one.fitness.values[1]))

    ax1.scatter(MA.MA_activator, list(ind_one['act']), alpha=0.5, c='b')
    ax2.scatter(MA.MA_inhibitor, list(ind_one['inh']), alpha=0.5, c='b')

    if ind_two is not None:
        ax1.scatter(MA.MA_activator, list(ind_two['act']), alpha=0.5, c='r')
        ax2.scatter(MA.MA_inhibitor, list(ind_two['inh']), alpha=0.5, c='r')
        ax1.set_title('Blue Score: '+str(ind_one.fitness.values[0])+'\n Red Score: '+str(ind_two.fitness.values[0]), fontsize=8)
        ax2.set_title('Blue Score: '+str(ind_one.fitness.values[1])+'\n Red Score: '+str(ind_two.fitness.values[1]), fontsize=8)

    return fig, (ax1, ax2)