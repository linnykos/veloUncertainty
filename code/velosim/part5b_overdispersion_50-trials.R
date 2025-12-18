rm(list=ls())
library(Matrix)

data_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/"

for(seed in 0:49){
  print(paste0("Working on seed: ", seed))
  
  S <- Matrix::readMM(paste0(data_folder,'spliced_seed', seed, '.mtx'))
  U <- Matrix::readMM(paste0(data_folder,'unspliced_seed', seed, '.mtx'))
  
  overdisp_S <- sapply(1:ncol(S), function(i) {
    res <- glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)
    1/res$overdispersions
  })
  overdisp_U <- sapply(1:ncol(U), function(i) {
    res <- glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)
    1/res$overdispersions
  })
  
  quantile(overdisp_S)
  quantile(overdisp_U)
  
  write.csv(overdisp_S, file=paste0(data_folder,'overdisp_S_seed', seed, '.csv'))
  write.csv(overdisp_U, file=paste0(data_folder,'overdisp_U_seed', seed, '.csv'))
}
