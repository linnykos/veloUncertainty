import numpy as np
from scipy.stats import gamma, beta
from scipy.sparse import csr_matrix, find
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf



# may set a global seed through np.random.seed()
################################################################
## This function is based on code from the Rcpp Gallery https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/.
## It draws one realization from a multinomial distribution with size parameter "size" and probability vector "probs".
## Returns a single multinomial realization, whose length is equal to the number of elements in "probs".
## @param size the multinomial size parameter. Should be a positive integer.
## @param probs the multinomial probability parameter. Should be a vector with non-negative entries that sum to 1.
def rmultinom_1(size, probs):
    """
    Draws one realization from a multinomial distribution with the specified size and probability vector.
    Args:
    size (int): The multinomial size parameter. Should be a positive integer.
    probs (list or np.array): The multinomial probability parameter. Should be a vector with non-negative entries that sum to 1.
    Returns:
    np.array: A single multinomial realization, whose length is equal to the number of elements in "probs".
    """
    if not isinstance(size, int) or size <= 0:
        raise ValueError("Size parameter should be a positive integer.")
    if not isinstance(probs, (list, np.ndarray)):
        raise ValueError("Probs should be a list or numpy array.")
    if not all(p >= 0 for p in probs):
        raise ValueError("All elements in probs should be non-negative.")
    if not np.isclose(sum(probs), 1):
        raise ValueError("The elements in probs should sum to 1.")
    outcome = np.random.multinomial(size, probs)
    return outcome

################################################################
## Returns a single Dirichlet-Multinomial realization, who length is equal to "folds".
## Assumes that we want all folds to have equal information in them- slower version must be called if this is not the case.
## The size parameter is "size", and the other parameters are all given by "b/folds".
## Corresponds to drawing a vector called "probs" from a Dirichlet distribution with all parameters equal to b/folds,
## and then drawing from a multinomial(size, probs) distribution.
def dir_mul_sample(size, folds, b):
    """
    Returns a single Dirichlet-Multinomial realization with the specified size and number of folds.
    Args:
    size (int): The multinomial size parameter.
    folds (int): The number of folds.
    b (float): The parameter for the Dirichlet distribution.
    Returns:
    np.array: A single Dirichlet-Multinomial realization.
    """
    # Initialize probs to equal probabilities for each fold
    probs = np.full(folds, 1.0 / folds)
    # If b is not infinity, sample "probs" from a Dirichlet distribution
    if not np.isinf(b):
        gammas = gamma.rvs(1.0 / folds * b, size=folds) # np.random.gamma(1.0 / folds * b, 1.0,size=folds)
        sum_gammas = np.sum(gammas)
        if sum_gammas > 0:
            probs = gammas / sum_gammas
        else:
            # Handle the case where sum_gammas is 0
            idx = np.random.choice(folds)
            probs = np.zeros(folds)
            probs[idx] = 1
    # Draw from a multinomial distribution with the calculated probs
    result = np.random.multinomial(size, probs)
    return result


################################################################
## Returns a single Dirichlet-Multinomial realization, who length is equal to "folds".
## Does not assume that all folds have equal information, so is slower than "dir_mul_sample_cpp".
## The size parameter is "size", and the other parameters are given by the vector "epsilon*b".
## Corresponds to drawing a vector called "probs" from a Dirichlet distribution with parameter "epsilon*b",
## and then drawing from a multinomial(size, probs) distribution.
def dir_mul_slower(size, epsilon, b):
    """
    Returns a single Dirichlet-Multinomial realization with the specified size and epsilon vector.
    Args:
    size (int): The multinomial size parameter.
    epsilon (list or np.array): The parameter vector for the Dirichlet distribution.
    b (float): The parameter for the Dirichlet distribution.
    Returns:
    np.array: A single Dirichlet-Multinomial realization.
    """
    # If b is infinity, we can use the epsilon vector directly
    if np.isinf(b):
        result = np.random.multinomial(size, epsilon)
        return result
    # Otherwise, initialize probs with length of epsilon
    folds = len(epsilon)
    probs = np.zeros(folds)
    # Generate gamma samples for each fold
    gammas = np.zeros(folds)
    for j in range(folds):
        gammas[j] = gamma.rvs(epsilon[j] * b, scale=1.0)
    sum_gammas = np.sum(gammas)
    # Normalize to get the probs vector
    if sum_gammas > 0:
        probs = gammas / sum_gammas
    else:
        # Handle the case where sum_gammas is 0
        idx = np.random.choice(folds, p=epsilon/np.sum(epsilon))
        probs[idx] = 1
    # Draw from a multinomial distribution with the calculated probs
    result = np.random.multinomial(size, probs)
    return result

################################################################
## Returns a single draw from a beta binomial distribution.
## This is the univariate version of the Dirichlet Multinomial,
## but implementing it separately allows for additional efficiency.
## Returns a single draw from BetaBinomial(x, epsilon*b, (1-epsilon)*b)
def beta_bin_sample(size, eps, b):
    """
    Returns a single draw from a Beta-Binomial distribution.
    Args:
    size (int): The binomial size parameter.
    eps (float): The epsilon parameter for the Beta distribution.
    b (float): The parameter for the Beta distribution.
    Returns:
    int: A single draw from the Beta-Binomial distribution.
    """
    # Special case: b is infinity
    if np.isinf(b):
        return np.random.binomial(size, eps)
    # Special case: b is 0
    if b == 0:
        choices = [size, 0]
        probs = [eps, 1 - eps]
        return np.random.choice(choices, p=probs)
    # General case
    p = beta.rvs(eps * b, (1.0 - eps) * b)
    return np.random.binomial(size, p)

################################################################
## Handles the efficient calling of dir_mul_sample across all elements in x
## Entry j in x gets turned into a vector with length folds, where these entries
## were drawn from a Dirichlet-multinomial distribution with parameters x and overdisps[j]/folds.
## Assumes that we are allocating information equally between folds.
## If this is not the case, the slower version must be called.
def mapply_dir_mul_sample(x, folds, overdisps):
    """
    Applies dir_mul_sample across all elements in x.
    Args:
    x (np.array or list): The integer vector of sizes.
    folds (int): The number of folds.
    overdisps (np.array or list): The overdispersion parameters.
    Returns:
    np.array: A matrix where each column is a Dirichlet-Multinomial realization.
    """
    n = len(x)
    result = np.zeros((folds, n), dtype=int)
    for i in range(n):
        sample = dir_mul_sample(x[i], folds, overdisps[i])
        result[:, i] = sample
    return result

################################################################
## Handles the efficient calling of dir_mul_sample across all elements in x.
## The length of the vector epsilon tells us how many folds we are generating.
## Entry j in x gets turned into a vector with length folds, where these entries
## were drawn from a Dirichlet-multinomial distribution with parameters x and epsilon*overdisps[j].
def mapply_dir_mul_slower(x, epsilon, overdisps):
    """
    Applies dir_mul_slower across all elements in x.
    Args:
    x (np.array or list): The integer vector of sizes.
    epsilon (np.array or list): The parameter vector for the Dirichlet distribution.
    overdisps (np.array or list): The overdispersion parameters.
    Returns:
    np.array: A matrix where each column is a Dirichlet-Multinomial realization.
    """
    n = len(x)
    folds = len(epsilon)
    result = np.zeros((folds, n), dtype=int)
    for i in range(n):
        sample = dir_mul_slower(x[i], epsilon, overdisps[i])
        result[:, i] = sample
    return result

################################################################
## Handles the efficient calling of beta_bin_sample_cpp across all elements in x.
## Entry j in x gets turned into (x1,x2), where x1 is drawn from a betabinomial distribution with parameters
## x, eps1*overdisp[j], and (1-eps)*overdisps[j]. And x2=x-x1.
def mapply_betabin_sample(x, eps1, overdisps):
    """
    Applies beta_bin_sample across all elements in x.
    Args:
    x (np.array or list): The integer vector of sizes.
    eps1 (float): The epsilon parameter for the Beta distribution.
    overdisps (np.array or list): The overdispersion parameters.
    Returns:
    np.array: A matrix where each column is a Beta-Binomial realization (x1, x - x1).
    """
    n = len(x)
    result = np.zeros((2, n), dtype=int)
    for i in range(n):
        x1 = beta_bin_sample(x[i], eps1, overdisps[i])
        result[0, i] = x1
        result[1, i] = x[i] - x1
    return result

################################################################
def countsplit(X, folds=2, epsilon=None, overdisps=None):
    if epsilon is None:
        epsilon = np.full(folds, 1.0 / folds)
    if overdisps is None:
        overdisps = np.full(X.shape[1], np.inf)
        print("As no overdispersion parameters were provided, Poisson count splitting will be performed.")
    if len(overdisps) != X.shape[1]:
        raise ValueError("You should provide one overdispersion parameter for every column of X. Make sure that your matrix X is cell-by-gene rather than gene-by-cell.")
    if len(epsilon) != folds or not np.isclose(np.sum(epsilon), 1) or np.any(epsilon <= 0):
        raise ValueError("The parameter epsilon should be a vector with length folds with positive entries that sum to 1.")
    if not isinstance(X, csr_matrix):
        X = csr_matrix(X)
    row, col, data = find(X)
    mapped_overdisps = overdisps[col]
    if not np.array_equal(epsilon, np.full(folds, 1.0 / folds)):
        if folds == 2:
            results = mapply_betabin_sample(data, epsilon[0], mapped_overdisps)
        else:
            results = mapply_dir_mul_slower(data, epsilon, mapped_overdisps)
    else:
        results = mapply_dir_mul_sample(data, folds, mapped_overdisps)
    partition = []
    for f in range(folds):
        Xfold = X.copy()
        Xfold.data = results[f, :].astype(float)
        partition.append(Xfold)
    return partition

################################################################
"""
Estimate overdispersion parameter
Input: a cell-by-gene matrix
Output: a vector of overdispersion parameters of length Ngenes
"""
def estimate_overdisps(X):
    res = []
    p = None
    if X.ndim==1: 
        p = 1
        data = { 'counts': X }
        df = pd.DataFrame(data)
        model = smf.negativebinomial('counts ~ 1', data=df)
        result = model.fit()
        res.append(result.params['alpha'])
    if X.ndim==2:
        p = X.shape[1]
        for col in range(p):
            y = X[:, col]
            data = { 'counts': y }
            df = pd.DataFrame(data)
            model = smf.negativebinomial('counts ~ 1', data=df)
            result = model.fit()
            res.append(result.params['alpha'])
    return res


