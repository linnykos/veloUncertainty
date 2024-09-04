library(Matrix)
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

S <- readMM(paste0(data_folder,'v2_pancreas/pan_spliced.mtx'))
U <- readMM(paste0(data_folder,'v2_pancreas/pan_unspliced.mtx'))

overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_folder,'v2_pancreas/pan_overdisp_S.csv'))
write.csv(overdisp_U, file=paste0(data_folder,'v2_pancreas/pan_overdisp_U.csv'))
