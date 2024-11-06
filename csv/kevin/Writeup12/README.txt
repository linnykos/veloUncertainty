pancreas_gene_likelihood.csv made by Yuhong on October 10, 2024 (and then corrected on Oct 16). These results are for the pancreas for scVelo.
========[vv Oct 10 vv]======
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
========[vv Oct 16 vv]======
velocity_genes' has both True and False
9:25
I guess it might be an intermediate variable in the gene filtering process
9:29
They are probably using these boolean values to filter out genes
9:29
In the filter_genes_dispersion function https://github.com/theislab/scvelo/blob/main/scvelo/preprocessing/utils.py#L228
utils.py
def filter_genes_dispersion(
<https://github.com/theislab/scvelo|theislab/scvelo>theislab/scvelo | Added by GitHub
9:30
I saw something like this
9:30
Since we followed the default subset=True, all the non-highly-variable genes are dropped
9:32
But this is different from the highly_variable_genes from the original adata object.
======


pancreas_utv_gene_fitloss.csv and pancreas_velovi_gene_velocityR2.csv were provided on November 4.
Yuhong's message:
===v====
For unitvelo: fit_loss, highly_variable, velocity_genes; for velovi (standard procedure with preprocess): velocity_r2, highly_variable, velocity_genes
6:15
VeloVI without preprocessing does not contain velocity_r2, since this is computed in the preprocess function.
6:16
'highly_variable' here are the computed values and are not from the raw dataset.
===^====