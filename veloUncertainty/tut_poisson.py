import numpy as np
from sklearn.cluster import KMeans
from statsmodels.api import GLM, families
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import scipy.stats as stats

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/veloUncertainty/tut/"

# https://anna-neufeld.github.io/countsplit.tutorials/articles/countsplit_tutorial.html
# simulate Poisson data with no true signal
np.random.seed(1)
n = 1000
p = 200
lambda_param = 5
X = np.random.poisson(lambda_param, size=(n, p)) # X.shape = (1000,200)

# The naive method
### First, we will demonstrate that estimating the clusters and testing for differential expression 
### using the same data does not control the Type 1 error rate. 
### We refer to the practice of using the same data for clustering and differential expression testing 
### as the “naive method” or “double dipping”.
# Step 1: Apply log transformation and perform k-means clustering
seed_kmeans = 1
log_X = np.log(X + 1)
kmeans = KMeans(n_clusters=2, random_state=seed_kmeans).fit(log_X)
clusters_full = kmeans.labels_  # len = 1000
# Step 2: Fit Poisson regression models and extract coefficients & pvalues
results_naive_coef = []
results_naive_pval = []
for i in range(X.shape[1]):
    u = X[:, i].astype(float)
    cluster_factor = pd.Categorical(clusters_full)
    design_matrix = pd.get_dummies(cluster_factor, drop_first=True)
    design_matrix = sm.add_constant(design_matrix)
    design_matrix = design_matrix.map(lambda x: 1 if x is True else (0 if x is False else x))
    design_matrix_array = design_matrix.values
    model = GLM(u, design_matrix_array, family=families.Poisson()).fit()
    coef = model.params[1]  # Extract the coefficient for the second cluster
    pval = model.pvalues[1]
    results_naive_coef.append(coef)
    results_naive_pval.append(pval)

results_naive_coef = np.array(results_naive_coef) # len = 200
results_naive_pval = np.array(results_naive_pval) # len = 200
np.quantile(results_naive_pval, [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])



# Create a Q-Q plot comparing sample quantiles to a uniform distribution
plt.clf()
fig, ax = plt.subplots()
stats.probplot(results_naive_pval, dist="uniform", plot=ax)
ax.get_lines()[1].set_color('red')  # Change the color of the reference line to red
ax.set_title('Q-Q Plot - Naive Method')
ax.set_xlabel('Theoretical Quantiles')
ax.set_ylabel('Sample Quantiles')
plt.savefig(fig_folder+"poisson_pval_naive.png")
plt.clf()


# Poisson countsplitting
np.random.seed(1)
n = 1000
p = 200
lambda_param = 5
X = np.random.poisson(lambda_param, size=(n, p)) # X.shape = (1000,200)
# .astype(float)

np.random.seed(1)
split = countsplit(X, folds=2, epsilon=np.array([0.5,0.5]))
Xtrain = split[0].toarray()
Xtest = split[1].toarray()
log_Xtrain = np.log(Xtrain + 1)
kmeans_countsplit = KMeans(n_clusters=2, random_state=seed_kmeans).fit(log_Xtrain)
clusters_train = kmeans_countsplit.labels_ 

results_countsplit_coef = []
results_countsplit_pval = []
for i in range(X.shape[1]):
    u = Xtest[:, i]
    cluster_train_factor = pd.Categorical(clusters_train)
    design_matrix = pd.get_dummies(cluster_train_factor, drop_first=True)
    design_matrix = sm.add_constant(design_matrix)
    design_matrix = design_matrix.map(lambda x: 1 if x is True else (0 if x is False else x))
    design_matrix_array = design_matrix.values
    model = GLM(u, design_matrix_array, family=families.Poisson()).fit()
    coef = model.params[1]
    pval = model.pvalues[1]
    results_countsplit_coef.append(coef)
    results_countsplit_pval.append(pval)

results_countsplit_coef = np.array(results_countsplit_coef) # len = 200
results_countsplit_pval = np.array(results_countsplit_pval) # len = 200
np.quantile(results_countsplit_pval, [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])

# Create a Q-Q plot comparing sample quantiles to a uniform distribution
plt.clf()
fig, ax = plt.subplots()
stats.probplot(results_countsplit_pval, dist="uniform", plot=ax)
ax.get_lines()[1].set_color('red')  # Change the color of the reference line to red
ax.set_title('Q-Q Plot - Countsplitting')
ax.set_xlabel('Theoretical Quantiles')
ax.set_ylabel('Sample Quantiles')
plt.savefig(fig_folder+"poisson_pval_countsplit.png")
plt.clf()

