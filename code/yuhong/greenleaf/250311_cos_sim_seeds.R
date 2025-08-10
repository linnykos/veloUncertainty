cos_sim.glf_GPC <- read.csv('/Users/y2564li/Downloads/proj_scRNA/data/cos_sim_5seeds_glf_GPC.csv')
cos_sim.glf_nGPCgrid <- read.csv('/Users/y2564li/Downloads/proj_scRNA/data/cos_sim_5seeds_glf_nGPCgrid.csv')

method = 'scv'

# median per seed
## scv
apply(cos_sim.glf_GPC[,grep('scv',names(cos_sim.glf_GPC))], MARGIN=2, FUN=median) %>% mean
apply(cos_sim.glf_GPC[,grep('scv',names(cos_sim.glf_GPC))], MARGIN=2, FUN=median) %>% sd
apply(cos_sim.glf_nGPCgrid[,grep('scv',names(cos_sim.glf_nGPCgrid))], MARGIN=2, FUN=median) %>% mean
apply(cos_sim.glf_nGPCgrid[,grep('scv',names(cos_sim.glf_nGPCgrid))], MARGIN=2, FUN=median) %>% sd

## utv
apply(cos_sim.glf_GPC[,grep('utv',names(cos_sim.glf_GPC))], MARGIN=2, FUN=median) %>% mean
apply(cos_sim.glf_GPC[,grep('utv',names(cos_sim.glf_GPC))], MARGIN=2, FUN=median) %>% sd
apply(cos_sim.glf_nGPCgrid[,grep('utv',names(cos_sim.glf_nGPCgrid))], MARGIN=2, FUN=median) %>% mean
apply(cos_sim.glf_nGPCgrid[,grep('utv',names(cos_sim.glf_nGPCgrid))], MARGIN=2, FUN=median) %>% sd


# mean per cell
apply(cos_sim.glf_GPC[,grep('scv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=mean) %>% hist(xlim=c(-1,1))
apply(cos_sim.glf_GPC[,grep('scv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=sd) %>% hist # sd

apply(cos_sim.glf_nGPCgrid[,grep('scv',names(cos_sim.glf_nGPCgrid))], MARGIN=1, FUN=mean) %>% hist(xlim=c(-1,1))
(apply(cos_sim.glf_GPC[,grep('scv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=mean) - 
    apply(cos_sim.glf_nGPCgrid[,grep('scv',names(cos_sim.glf_nGPCgrid))], MARGIN=1, FUN=mean)) %>% hist(main='GPC - nGPC')

apply(cos_sim.glf_GPC[,grep('utv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=mean) %>% hist(xlim=c(-1,1))
apply(cos_sim.glf_nGPCgrid[,grep('utv',names(cos_sim.glf_nGPCgrid))], MARGIN=1, FUN=mean) %>% hist(xlim=c(-1,1))

(apply(cos_sim.glf_GPC[,grep('utv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=mean) - 
  apply(cos_sim.glf_nGPCgrid[,grep('utv',names(cos_sim.glf_nGPCgrid))], MARGIN=1, FUN=mean)) %>% hist(main='GPC - nGPC')

apply(cos_sim.glf_GPC[,grep('utv',names(cos_sim.glf_GPC))], MARGIN=1, FUN=sd) %>% hist # sd

############

pvalue_vec <- sapply(2:6, function(j){
  t.test(cos_sim.glf_GPC[,j],cos_sim.glf_nGPCgrid[,j], paired = TRUE, alternative = "greater")$p.value
})

pvalue_vec <- sapply(7:11, function(j){
  t.test(cos_sim.glf_GPC[,j],cos_sim.glf_nGPCgrid[,j], paired = TRUE, alternative = "greater")$p.value
})

plot(cos_sim.glf_GPC[,9], cos_sim.glf_nGPCgrid[,9], asp = TRUE, pch = 16, col = rgb(0.5,0.5,0.5,0.5))

plot(cos_sim.glf_GPC[,2], cos_sim.glf_nGPCgrid[,2], asp = TRUE, pch = 16, col = rgb(0.5,0.5,0.5,0.01))
