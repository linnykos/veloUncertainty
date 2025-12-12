rm(list=ls())
library(Matrix)

data_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/"

S <- Matrix::readMM(paste0(data_folder,'spliced.mtx'))
U <- Matrix::readMM(paste0(data_folder,'unspliced.mtx'))

overdisp_S <- sapply(1:ncol(S), function(i) {
  res <- glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)
  1/res$overdispersions
})
overdisp_U <- sapply(1:ncol(U), function(i) {
  res <- glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)
  1/res$overdispersions
})

write.csv(overdisp_S, file=paste0(data_folder,'overdisp_S.csv'))
write.csv(overdisp_U, file=paste0(data_folder,'overdisp_U.csv'))