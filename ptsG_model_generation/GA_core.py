# Modules
import numpy as np
import pandas as pd

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
                   n_iter = int,
                   grid = pd.Series,
                   stats = None,
                   hall_of_fame = None,
                   verbose = __debug__):
    """
    Modified DEAP mu+lambda evolutionary algorithm using var_or

    :param list pop: List of individuals to serve as the starting pop
    :param object rng: np.random.generator or np.random.RandomState class
    :param base.Toolbox() toolbox: DEAP class containing evolution operators
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