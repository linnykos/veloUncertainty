import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
from scipy.stats import nbinom
from sklearn.cluster import KMeans
from statsmodels.api import GLM, families
import matplotlib.pyplot as plt
import scipy.stats as stats


import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/veloUncertainty/tut/"

np.random.seed(1)
n = 1000  
p = 200   
mu = 5    # Mean of the negative binomial distribution
size = 5  # Size parameter of the negative binomial distribution
# Generate the negative binomial distributed data
data = nbinom.rvs(n=size, p=size/(size + mu), size=n*p)
# Reshape the data into a matrix with n rows and p columns
X = data.reshape((n, p))

seed_kmeans=1

def get_pvalues(Xtrain,Xtest,seed_kmeans=1):
    kmeans = KMeans(n_clusters=2, random_state=seed_kmeans).fit(np.log(Xtrain+1))
    clusters_full = kmeans.labels_ 
    results_pval = []
    for i in range(X.shape[1]):
        u = Xtest[:, i]
        cluster_factor = pd.Categorical(clusters_full)
        design_matrix = pd.get_dummies(cluster_factor, drop_first=True)
        design_matrix = sm.add_constant(design_matrix)
        design_matrix = design_matrix.map(lambda x: 1 if x is True else (0 if x is False else x))
        design_matrix_array = design_matrix.values
        model = GLM(u, design_matrix_array, family=families.Poisson()).fit()
        pval = model.pvalues[1]
        results_pval.append(pval)
    return results_pval

def plot_pvalues(Xtrain, Xtest, method, fig_name, seed_kmeans=1):
    pval = get_pvalues(Xtrain, Xtest, seed_kmeans)
    plt.clf()
    fig, ax = plt.subplots()
    stats.probplot(pval, dist="uniform", plot=ax)
    ax.get_lines()[1].set_color('red')  # Change the color of the reference line to red
    ax.set_title('Q-Q Plot - '+method)
    ax.set_xlabel('Theoretical Quantiles')
    ax.set_ylabel('Sample Quantiles')
    plt.savefig(fig_folder+fig_name)
    plt.clf()

get_pvalues(X,X)
plot_pvalues(X,X,method="naive",fig_name="NB_naive.png")
np.quantile(get_pvalues(X,X), [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])

## Poisson countsplitting
np.random.seed(1)
split = countsplit(X)
Xtrain = split[0].toarray()
Xtest = split[1].toarray()
get_pvalues(Xtrain,Xtest)
plot_pvalues(Xtrain,Xtest,method="countsplit",fig_name="NB_countsplitPoi.png")

## NB countsplitting
np.random.seed(1)
splitNB = countsplit(X, overdisps=np.full(p,size)) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
#get_pvalues(XtrainNB,XtestNB)
#np.quantile(get_pvalues(XtrainNB,XtestNB), [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB)",fig_name="NB_countsplitNB.png")

## NB countsplitting with correct overdispersion estimation
np.random.seed(1)
overdisps = estimate_overdisps(X)
splitNB = countsplit(X, overdisps=overdisps) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB size=5, overdisps estimated)",fig_name="NB_countsplitNB_size5_overdisps.png")

## NB countsplitting with wrong overdispersion estimation
np.random.seed(1)
overdisps = estimate_overdisps(X)
splitNB = countsplit(X, overdisps=1/overdisps) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB size=5, overdisps estimated but use the inverse)",fig_name="NB_countsplitNB_size5_overdispsInverse.png")

### size=.05
np.random.seed(1)
n = 1000  
p = 200   
mu = 5    # Mean of the negative binomial distribution
size = .05  # Size parameter of the negative binomial distribution
# Generate the negative binomial distributed data
data = nbinom.rvs(n=size, p=size/(size + mu), size=n*p)
# Reshape the data into a matrix with n rows and p columns
X = data.reshape((n, p))

plot_pvalues(X,X,method="naive",fig_name="NB_size5e-2_naive.png")

seed_kmeans=1
## Poisson countsplitting
np.random.seed(1)
split = countsplit(X)
Xtrain = split[0].toarray()
Xtest = split[1].toarray()
#get_pvalues(Xtrain,Xtest)
plot_pvalues(Xtrain,Xtest,method="countsplit",fig_name="NB_size5e-2_countsplitPoi.png")

## NB countsplitting
np.random.seed(1)
splitNB = countsplit(X, overdisps=np.full(p,size)) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
#get_pvalues(XtrainNB,XtestNB)
#np.quantile(get_pvalues(XtrainNB,XtestNB), [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB size=.05)",fig_name="NB_size5e-2_countsplitNB.png")

## NB countsplitting with correct overdispersion estimation
np.random.seed(1)
overdisps = estimate_overdisps(X)
splitNB = countsplit(X, overdisps=overdisps) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB size=.05, overdisps estimated)",fig_name="NB_size5e-2_countsplitNB_overdisps.png")

## NB countsplitting with wrong overdispersion estimation
np.random.seed(1)
overdisps = estimate_overdisps(X)
splitNB = countsplit(X, overdisps=1/overdisps) # size=1/alpha, where alpha is the obtained overdispersion parameter estimate from estimation_overdisps()
XtrainNB = splitNB[0].toarray()
XtestNB = splitNB[1].toarray()
plot_pvalues(XtrainNB,XtestNB,method="countsplit (NB size=.05, overdisps estimated but use the inverse)",fig_name="NB_size5e-2_countsplitNB_overdispsInverse.png")



######### testing code #########
####################################
np.random.seed(42)
n = 1000000  # number of observations
X = np.random.normal(size=n)
# True parameters
mu = [1] * int(n)
mu = np.full(len(mu), 1, dtype=int)
size = 1/5  # size is equal to 1/(over-dispersion) parameter
y = np.random.negative_binomial(size, size / (size + mu))

np.mean(y)
np.var(y)
r = size
p = size / (size+mu[0])
r * (1-p) / p # this should be equal to the true mean
r * (1-p) / p**2 # this should equal the true variance

# data = np.column_stack((y, mu))
# model = sm.GLM(y, data, family=sm.families.NegativeBinomial()).fit()

data = { 'counts': y }
df = pd.DataFrame(data)

model = smf.negativebinomial('counts ~ 1', data=df)
result = model.fit()
overdispersion_parameter = result.params['alpha']
print(overdispersion_parameter)

#########




