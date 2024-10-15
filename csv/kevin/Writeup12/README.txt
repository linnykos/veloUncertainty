pancreas_gene_likelihood.csv made by Yuhong on October 10. These results are for the pancreas for scVelo.
========
for each total/split, I include a column 'highly_variable' and a 'lik'. 'highly_variable' is copied from the var['highly_variable'] and 'lik' is for 'fit_likelihood'.
8:55
Column prefixes are named in the form of split seed + split#, e.g. 317s1 means seed=317, split1
8:57
If a gene is contained in a split, the 'highly_variable' column is either True or False. If it is not included, the column is NaN.
8:58
Some genes are included in a split but have no likelihood value. They can be detected by finding a 'highly_variable' value (T/F) but NaN in 'lik'.
8:59
Last 4 columns are overdispersion parameter estimates and nonzero%, separated for spliced and unspliced.
Changed high_variable columns into 1(True)/0(False)/NaN
========