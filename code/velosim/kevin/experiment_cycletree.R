rm(list=ls())

library(VeloSim)

# number of total cells
ncells_total <- 1000
# number of cell in cell-cycle stage
ncells_cycle <- 250
# number of total genes
ngenes <- 1000
# total number of evfs: homo-evf + diff-evf
nevf <- 20
# number of diff-evf
n_de_evf <- 12
# which kinetic parameters is changed according to the sampling position, n_de_evf is used in s
vary <- "all"
# the standard deviation of the gaussian distribution, denote the noise size
Sigma <- 0.4
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
randseed <- 4
# number of cells to be removed at the beginning of the simulation
n_unstable <- 10
# number of cells in the second cycle
nextra <- 0

# read in tree
phyla <- ape::read.tree(system.file("extdata", "Newick_ABCD.txt", package = "VeloSim"))
par(mfrow=c(1,2))
plot(phyla)

result_true <- SimulateCycleTree(ncells_total=ncells_total, ncells_cycle = ncells_cycle, ngenes=ngenes, evf_center=1,nevf=nevf,
                                 randseed=randseed, n_de_evf=n_de_evf,vary=vary,Sigma=Sigma, phyla = phyla,
                                 geffect_mean=geffect_mean,gene_effects_sd=gene_effects_sd,
                                 gene_effect_prob=gene_effect_prob,bimod=bimod,param_realdata=param_realdata,
                                 scale_s=scale_s, prop_hge=0.015, mean_hge=5, n_unstable=n_unstable, nextra=nextra)
result <- result_true

data_dir <- '/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/cycleTree/'
# write.csv(result$counts_u, file = paste0(data_dir,"velosim_counts_u.csv"), row.names=F)
# write.csv(result$counts_s, file = paste0(data_dir,"velosim_counts_s.csv"), row.names=F)
# write.csv(result$cell_time, file = paste0(data_dir,"velosim_pseudo_time.csv"), row.names=F)
# write.csv(result$velocity, file = paste0(data_dir,"velosim_velocity.csv"), row.names=F)

setwd(data_dir)
# Plot with cell developmental time(pseudo-time), default PCA
plotPseudotime(filename = "pseudotime.pdf", result, width = 6, height = 5)
# Plot with RNA velocity
plotVelo(filename = "velocity.pdf", result, arrow.length = 5, width = 6, height = 5)
# Plot with cell developmental time(pseudo-time), using UMAP
plotPseudotime(filename = "pseudotime_umap.pdf", result, dr = "umap", width = 6, height = 5)


