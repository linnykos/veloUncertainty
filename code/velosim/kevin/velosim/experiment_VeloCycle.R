rm(list=ls())

# https://github.com/PeterZZQ/VeloSim/blob/1a7cf83eb5b032a88a020cb060ae1f76b98fcd01/DESCRIPTION
library(plyr)
library(ggplot2)
library(umap)

source("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/code/velosim/kevin/custom_simulate_kinetics.R")
source("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/code/velosim/kevin/custom_simulate_counts.R")
source("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/code/velosim/kevin/custom_plotting.R")

# 1. siulation ####
ncells_total     <- 1000
ngenes           <- 1000

nevf             <- 20
n_de_evf         <- 12
vary             <- "all"      # fine for now
Sigma            <- 0.4
evf_center       <- 1

gene_effect_prob <- 0.3
geffect_mean     <- 0
gene_effects_sd  <- 1

param_realdata   <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/param_realdata.zeisel.imputed_KZL.RData"
bimod            <- 0

scale_s          <- 1
prop_hge         <- 0.015
mean_hge         <- 5

nextra           <- 100
n_unstable       <- 10
randseed         <- 3

## rate=1 ####
set.seed(1126)
# simulate data
result <- SimulateVeloCycle(
  ncells_total = ncells_total,
  ngenes       = ngenes,
  evf_center   = evf_center,
  nevf         = nevf,
  randseed     = randseed,
  n_de_evf     = n_de_evf,
  vary         = vary,
  Sigma        = Sigma,
  geffect_mean = geffect_mean,
  gene_effects_sd   = gene_effects_sd,
  gene_effect_prob  = gene_effect_prob,
  bimod        = bimod,
  param_realdata    = param_realdata,
  scale_s      = scale_s,
  prop_hge     = prop_hge,
  mean_hge     = mean_hge,
  nextra       = nextra,
  n_unstable   = n_unstable,
  plot         = FALSE
)

###########

data_dir <- '/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/'
write.csv(result$counts_u, file = paste0(data_dir,"velosim_counts_u.csv"), row.names=F)
write.csv(result$counts_s, file = paste0(data_dir,"velosim_counts_s.csv"), row.names=F)
write.csv(result$cell_time, file = paste0(data_dir,"velosim_pseudo_time.csv"), row.names=F)
write.csv(result$velocity, file = paste0(data_dir,"velosim_velocity.csv"), row.names=F)

cur_dir <- getwd()
setwd(data_dir)
# Plot with cell developmental time(pseudo-time), default PCA
plotPseudotime(filename = "pseudotime.pdf", result, width = 6, height = 5)
# Plot with RNA velocity
plotVelo(filename = "velocity.pdf", result, width = 6, height = 5,
         arrow.length = 2)
# Plot with cell developmental time(pseudo-time), using UMAP
plotPseudotime(filename = "pseudotime_umap.pdf", result, dr = "umap", width = 6, height = 5)
setwd(cur_dir)
