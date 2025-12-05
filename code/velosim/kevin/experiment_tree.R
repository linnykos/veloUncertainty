rm(list=ls())

library(VeloSim)

# number of total cells
ncells_total <- 2000
# number of total genes
ngenes <- 500
# total number of evfs: homo-evf + diff-evf
nevf <- 20
# number of diff-evf
n_de_evf <- 12
# which kinetic parameters is changed according to the sampling position, n_de_evf is used in s
vary <- "all"
# the standard deviation of the gaussian distribution, denote the noise size
Sigma <- 0.1
# the mean of the homo-evf, gaussian distribution again
evf_center <- 1
# with which probability the gene effect value is preserved
gene_effect_prob <- 0.3
# the mean of the gaussian distribution of gene effect
geffect_mean <- 0
# the standard deviation of the gaussian distribution of gene effect
gene_effects_sd <- 1
# the real data to obtain the kinetic parameter distribution from
param_realdata <- "zeisel.imputed"

bimod <- 0
scale_s <- 1
randseed <- 2

# read in tree
phyla <- ape::read.tree(system.file("extdata", "Newick_ABCD.txt", package = "VeloSim"))
par(mfrow=c(1,2))
plot(phyla)

result <- SimulateVeloTree(ncells_total=ncells_total,ngenes=ngenes, evf_center=1,nevf=nevf,
                           phyla=phyla, randseed=randseed, n_de_evf=n_de_evf,vary=vary,Sigma=Sigma,geffect_mean=geffect_mean,
                           gene_effects_sd=gene_effects_sd,gene_effect_prob=gene_effect_prob,
                           bimod=bimod,param_realdata=param_realdata,scale_s=scale_s,
                           prop_hge=0.015, mean_hge=5, n_unstable=0, plot = FALSE)

data_dir <- '/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/tree/'

# Plot with cell developmental time(pseudo-time), default PCA
cur_dir <- getwd()
setwd(data_dir)
plotPseudotime(file = "pseudotime.pdf", result, width = 6, height = 5)
# Plot with RNA velocity
plotVelo(filename = "velocity.pdf", result, arrow.length = 2, width = 6, height = 5)
# Plot with cell developmental time(pseudo-time), using UMAP
plotPseudotime(filename = "pseudotime_umap.pdf", result, dr = "umap", width = 6, height = 5)
setwd(cur_dir)
