library(glmGamPoi)

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

dim(S)
dim(U)
fit <- glm_gp(S[,1], design = ~ 1)$overdispersions

apply(as.matrix(S[,11]),MARGIN=2, function(x) {glmGamPoi::glm_gp(dataframe(x=as.numeric(x)), design=~1)$overdispersions})

glmGamPoi::glm_gp(as.matrix(x=as.numeric(S[,11])), design=~1)

apply(S,MARGIN=2, function(x) {glmGamPoi::glm_gp(x, design=~1)$overdispersions}) # MARGIN=2 -> columns


library(Seurat)
library(SeuratDisk)
library(countsplit)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

scvelo_sc <- zellkonverter::readH5AD(file=paste0(data_folder,'v2_larry/larry_total_allgenes.h5ad'))

S <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced", layer = "data")
U <- SeuratObject::LayerData(scvelo_seurat, assay = "unspliced", layer = "data")
dim(S)
dim(U)

library(Matrix)

# Load the Matrix Market file
sparse_matrix <- readMM(paste0(data_folder,'v2_larry/S_tmp.mtx'))
S = sparse_matrix

fit <- glm_gp(S[,1], design = ~ 1)$overdispersions
apply(S, MARGIN=2, function(u) { MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(u)))$theta })

df = data.frame(x = as.numeric(S[,1]))
MASS::glm.nb(x ~ 1, data=df)$theta
# [6] err, err, 0.0001359058, 0.009956343, 0.9398475
# [11] 0.6428239, err, 0.0840699, 0.8152222, 0.9489258

x = S[,5]
df = data.frame(x = as.numeric(x))
MASS::glm.nb(x ~ 1, data=df)$theta

overdisp = MASS::glm.nb(x ~ 1, data=df)$theta
mean(x)
var(x)
mean(x)+mean(x)^2/overdisp
mean(x)+mean(x)^2*overdisp

apply(S[,1:15],MARGIN=2, function(x) {glmGamPoi::glm_gp(x, design=~1)$overdispersions})

apply(S,MARGIN=2, function(x) {glmGamPoi::glm_gp(x, design=~1)$overdispersions}) # MARGIN=2 -> columns

# S[,1:5]: 3.273167e-05, error, error, 0.536813, error
#2: Error in while ((it <- it + 1) < limit && abs(del) > eps) { : 
#  missing value where TRUE/FALSE needed
#In addition: Warning message:
#glm.fit: algorithm did not converge 

options(matrixStats.useNames = FALSE)
counts <- rnbinom(n = 10, mu = 5, size = 1/0.7)
fit <- glm_gp(counts, design = ~ 1)
fit
as.list(fit)[1:2]


set.seed(2410)
samples <- matrix(rpois(1000, lambda = 5),nrow=200,ncol=5)
apply(samples, MARGIN=2, function(u) { MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(u)))$theta })
# 54.20978 40543.89042 63599.66437 62264.73390 26969.01353

##########
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
df = read.csv(paste0(data_folder,'v2_larry/larry_df_allgenes.csv'))






