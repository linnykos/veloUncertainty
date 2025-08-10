library(glmGamPoi)
library(Matrix)

#from scipy.io import mmwrite
#mmwrite('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/greenleaf_spliced.mtx', adata.layers['spliced'])
#mmwrite('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/greenleaf_unspliced.mtx', adata.layers['unspliced'])


## estimate overdispersion parameter using glm_gp

# Load the sparse matrix from the Matrix Market file
S <- readMM("/Users/y2564li/Downloads/greenleaf_spliced.mtx")
U <- readMM("/Users/y2564li/Downloads/greenleaf_unspliced.mtx")

# comparison
MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(S[,1])))$theta
MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(S[,3])))$theta
glmGamPoi::glm_gp(S[,2], design = ~ 1)$overdispersions


overdisp_S <- sapply(1:ncol(S), function(i) 1/glmGamPoi::glm_gp(as.numeric(S[,i]), design = ~ 1)$overdispersions)

overdisp_U <- sapply(1:ncol(U), function(i) 1/glmGamPoi::glm_gp(as.numeric(U[,i]), design = ~ 1)$overdispersions)

"
overdisp_S <- apply(S, MARGIN=2, function(u) {
  print('*')
  if (sum(u)==0) Inf
  else 1/glmGamPoi::glm_gp(as.numeric(u), design = ~ 1)$overdispersions })

overdisp_U <- apply(U, MARGIN=2, function(u) {
  print('*')
  if (sum(u)==0) Inf
  else 1/glmGamPoi::glm_gp(as.numeric(u), design = ~ 1)$overdispersions })
"

write.csv(overdisp_S, file='/Users/y2564li/Downloads/proj_scRNA/greenleaf_overdisp_S.csv')
write.csv(overdisp_U, file='/Users/y2564li/Downloads/proj_scRNA/greenleaf_overdisp_U.csv')


# Note from larry: 
## when there is only 1 count, e.g. the first gene, glmGamPoi returns 0 where the inverse is Inf, while nb gives 3.273167e-05


