library(Matrix)
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

dataset_short = 'larryMult'

S <- readMM(paste0(data_folder,'v4_',dataset_short,'/',dataset_short,'_spliced.mtx'))
U <- readMM(paste0(data_folder,'v4_',dataset_short,'/',dataset_short,'_unspliced.mtx'))

overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_folder,'v4_',dataset_short,'/',dataset_short,'_overdisp_S.csv'))
write.csv(overdisp_U, file=paste0(data_folder,'v4_',dataset_short,'/',dataset_short,'_overdisp_U.csv'))
