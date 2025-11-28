library(Matrix)
library(glmGamPoi)
data_folder = "/Users/wakeup/Downloads/"

S = read.csv(paste0(data_folder,'velosim_counts_s.csv'))
U = read.csv(paste0(data_folder,'velosim_counts_u.csv'))

overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)
overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

write.csv(overdisp_S, file=paste0(data_folder,'velosim_overdisp_S.csv'), row.names=F)
write.csv(overdisp_U, file=paste0(data_folder,'velosim_overdisp_U.csv'), row.names=F)



