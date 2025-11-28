devtools::load_all()
library(VeloSim)

# 1. siulation ####
# number of total cells
ncells_total <- 500
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
# param_realdata <- NULL

bimod <- 0
scale_s <- 1
randseed <- 3
# number of cells in the second cycle
nextra <- 100
# number of cells to be removed at the beginning of the simulation
n_unstable <- 10

## rate=1 ####
set.seed(1126)
# simulate data
result_true <- SimulateVeloCycle(ncells_total=ncells_total,ngenes=ngenes, evf_center=1,nevf=nevf, 
                                 randseed=randseed, n_de_evf=n_de_evf,vary=vary,Sigma=Sigma,
                                 geffect_mean=geffect_mean,gene_effects_sd=gene_effects_sd,
                                 gene_effect_prob=gene_effect_prob,bimod=bimod,param_realdata=param_realdata,
                                 scale_s=scale_s, prop_hge=0.015, mean_hge=5, nextra = nextra, n_unstable=n_unstable, 
                                 plot = TRUE)
# add technical noise
result <- technicalNoise(result_true, capture.rate = 1)

data_dir <- '/Users/wakeup/Downloads/proj_scRNA/2511_simdata/251126-1_rate1/'
write.csv(result$counts_u, file = paste0(data_dir,"velosim_counts_u.csv"), row.names=F)
write.csv(result$counts_s, file = paste0(data_dir,"velosim_counts_s.csv"), row.names=F)
write.csv(result$cell_time, file = paste0(data_dir,"velosim_pseudo_time.csv"), row.names=F)
write.csv(result$velocity, file = paste0(data_dir,"velosim_velocity.csv"), row.names=F)

setwd(data_dir)
# Plot with cell developmental time(pseudo-time), default PCA
plotPseudotime(filename = "pseudotime.pdf", result, width = 6, height = 5)
# Plot with RNA velocity
plotVelo(filename = "velocity.pdf", result, arrow.length = 0.3, width = 6, height = 5)
# Plot with cell developmental time(pseudo-time), using UMAP
plotPseudotime(filename = "pseudotime_umap.pdf", result, dr = "umap", width = 6, height = 5)

## rate=0.8 ####
set.seed(1126)
# simulate data
result_true <- SimulateVeloCycle(ncells_total=ncells_total,ngenes=ngenes, evf_center=1,nevf=nevf, 
                                 randseed=randseed, n_de_evf=n_de_evf,vary=vary,Sigma=Sigma,
                                 geffect_mean=geffect_mean,gene_effects_sd=gene_effects_sd,
                                 gene_effect_prob=gene_effect_prob,bimod=bimod,param_realdata=param_realdata,
                                 scale_s=scale_s, prop_hge=0.015, mean_hge=5, nextra = nextra, n_unstable=n_unstable, 
                                 plot = TRUE)
# add technical noise
result <- technicalNoise(result_true, capture.rate = .8)

data_dir <- '/Users/wakeup/Downloads/proj_scRNA/2511_simdata/251126-2_rate.8/'
write.csv(result$counts_u, file = paste0(data_dir,"velosim_counts_u.csv"), row.names=F)
write.csv(result$counts_s, file = paste0(data_dir,"velosim_counts_s.csv"), row.names=F)
write.csv(result$cell_time, file = paste0(data_dir,"velosim_pseudo_time.csv"), row.names=F)
write.csv(result$velocity, file = paste0(data_dir,"velosim_velocity.csv"), row.names=F)

setwd(data_dir)
# Plot with cell developmental time(pseudo-time), default PCA
plotPseudotime(filename = "pseudotime.pdf", result, width = 6, height = 5)
# Plot with RNA velocity
plotVelo(filename = "velocity.pdf", result, arrow.length = 0.3, width = 6, height = 5)
# Plot with cell developmental time(pseudo-time), using UMAP
plotPseudotime(filename = "pseudotime_umap.pdf", result, dr = "umap", width = 6, height = 5)


# 2. overdispersion ####
library(Matrix)
library(glmGamPoi)
## rate=1 ####
data_dir <- '/Users/wakeup/Downloads/proj_scRNA/2511_simdata/251126-1_rate1/'
S = read.csv(paste0(data_dir,'velosim_counts_s.csv'))
U = read.csv(paste0(data_dir,'velosim_counts_u.csv'))

overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_dir,'velosim_overdisp_S.csv'), row.names=F)
write.csv(overdisp_U, file=paste0(data_dir,'velosim_overdisp_U.csv'), row.names=F)


## rate=0.8 ####
data_dir <- '/Users/wakeup/Downloads/proj_scRNA/2511_simdata/251126-2_rate.8/'
S = read.csv(paste0(data_dir,'velosim_counts_s.csv'))
U = read.csv(paste0(data_dir,'velosim_counts_u.csv'))

overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_dir,'velosim_overdisp_S.csv'), row.names=F)
write.csv(overdisp_U, file=paste0(data_dir,'velosim_overdisp_U.csv'), row.names=F)








