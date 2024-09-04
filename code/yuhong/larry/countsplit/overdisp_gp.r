library(Matrix)
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

larry_S <- readMM(paste0(data_folder,'v2_larry/larry_spliced.mtx'))
larry_U <- readMM(paste0(data_folder,'v2_larry/larry_unspliced.mtx'))

overdisp_S <- sapply(1:ncol(larry_S), function(i) 1/glmGamPoi::glm_gp(as.numeric(larry_S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(larry_U), function(i) 1/glmGamPoi::glm_gp(as.numeric(larry_U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_folder,'v2_larry/larry_overdisp_S.csv'))
write.csv(overdisp_U, file=paste0(data_folder,'v2_larry/larry_overdisp_U.csv'))

# test
# sapply(20:24, function(i) 1/glmGamPoi::glm_gp(as.numeric(larry_S[,i]), design = ~ 1)$overdispersions)