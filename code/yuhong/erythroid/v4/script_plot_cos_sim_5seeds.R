ery_Mark <- read.csv('/Users/y2564li/Downloads/proj_scRNA/pancreas/erythroid/cos_sim_5seeds_ery_Mark.csv')
ery_nMark <- read.csv('/Users/y2564li/Downloads/proj_scRNA/pancreas/erythroid/cos_sim_5seeds_ery_nMark.csv')
ery_nMark227 <- read.csv('/Users/y2564li/Downloads/proj_scRNA/pancreas/erythroid/cos_sim_5seeds_ery_nMark227.csv')

pan_Mark <- read.csv('/Users/y2564li/Downloads/proj_scRNA/pancreas/cos_sim_5seeds_pan_Mark.csv')
pan_nMark <- read.csv('/Users/y2564li/Downloads/proj_scRNA/pancreas/cos_sim_5seeds_pan_nMark.csv')

glf_GPC <- read.csv('/Users/y2564li/Downloads/proj_scRNA/greenleaf/cos_sim_5seeds_glf_GPC.csv')
glf_nGPCblk <- read.csv('/Users/y2564li/Downloads/proj_scRNA/greenleaf/cos_sim_5seeds_glf_nGPCblk.csv')

# erythorid plot
plot(ery_Mark[, names(ery_Mark)[grep('scv', names(ery_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='scv_Mark' )
plot(ery_Mark[, names(ery_Mark)[grep('utv', names(ery_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='utv_Mark' )
plot(ery_Mark[, names(ery_Mark)[grep('sct', names(ery_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='sct_Mark' )
plot(ery_Mark[, names(ery_Mark)[grep('velovi_seed', names(ery_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='velovi_Mark' )
plot(ery_Mark[, names(ery_Mark)[grep('velovi_woprep', names(ery_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='velovi_woprep_Mark' )

# mean and median for 5 seeds' cosine similarity
rbind( Mark_median=apply(ery_Mark[, names(ery_Mark)[grep('scv', names(ery_Mark))] ], 2, median),
       nMark_median=apply(ery_nMark[, names(ery_nMark)[grep('scv', names(ery_nMark))] ], 2, median),
       nMark227_median=apply(ery_nMark227[, names(ery_nMark227)[grep('scv', names(ery_nMark227))] ], 2, median))

rbind( Mark_mean=apply(ery_Mark[, names(ery_Mark)[grep('scv', names(ery_Mark))] ], 2, mean),
       nMark_mean=apply(ery_nMark[, names(ery_nMark)[grep('scv', names(ery_nMark))] ], 2, mean))




# pancreas
plot(pan_Mark[, names(pan_Mark)[grep('scv', names(pan_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='scv_Mark' )
plot(pan_Mark[, names(pan_Mark)[grep('utv', names(pan_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='utv_Mark' )
plot(pan_Mark[, names(pan_Mark)[grep('sct', names(pan_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='sct_Mark' )
plot(pan_Mark[, names(pan_Mark)[grep('velovi_seed', names(pan_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='velovi_Mark' )
plot(pan_Mark[, names(pan_Mark)[grep('velovi_woprep', names(pan_Mark))] ], col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2), 
     cex=.1, xlim=c(-1,1),ylim=c(-1,1), asp=T, main='velovi_woprep_Mark' )

rbind( Mark_median=apply(pan_Mark[, names(pan_Mark)[grep('scv', names(pan_Mark))] ], 2, median),
       nMark_median=apply(pan_nMark[, names(pan_nMark)[grep('scv', names(pan_nMark))] ], 2, median) )

apply(pan_Mark[, names(pan_Mark)[grep('scv', names(pan_Mark))] ], 2, median)
apply(pan_Mark[, names(pan_Mark)[grep('utv', names(pan_Mark))] ], 2, median)
apply(pan_Mark[, names(pan_Mark)[grep('sct', names(pan_Mark))] ], 2, median)
apply(pan_Mark[, names(pan_Mark)[grep('velovi_seed', names(pan_Mark))] ], 2, median)
apply(pan_Mark[, names(pan_Mark)[grep('velovi_woprep', names(pan_Mark))] ], 2, median)

apply(pan_nMark[, names(pan_nMark)[grep('scv', names(pan_nMark))] ], 2, median)
apply(pan_nMark[, names(pan_nMark)[grep('utv', names(pan_nMark))] ], 2, median)
apply(pan_nMark[, names(pan_nMark)[grep('sct', names(pan_nMark))] ], 2, median)
apply(pan_nMark[, names(pan_nMark)[grep('velovi_seed', names(pan_nMark))] ], 2, median)
apply(pan_nMark[, names(pan_nMark)[grep('velovi_woprep', names(pan_nMark))] ], 2, median)


# glf
rbind( GPC_median=apply(glf_GPC[, names(glf_GPC)[grep('scv', names(glf_GPC))] ], 2, median),
       nGPCblk_median=apply(glf_nGPCblk[, names(glf_nGPCblk)[grep('scv', names(glf_nGPCblk))] ], 2, median) )

###############

method_vec <- c("scv", "utv", "sct", "velovi_seed", "velovi_woprep")
teststat_vec <- rep(NA, length(method_vec))
names(teststat_vec) <- method_vec
for(i in 1:length(method_vec)){
  method_name <- method_vec[i]
  marker_values <- apply(pan_Mark[, names(pan_Mark)[grep(method_name, names(pan_Mark))] ], 2, median)
  nmarker_values <- apply(pan_nMark[, names(pan_nMark)[grep(method_name, names(pan_nMark))] ], 2, median)
  
  numerator <- mean(marker_values) - mean(nmarker_values)
  denominator <- sqrt(var(marker_values)/length(marker_values) + var(nmarker_values)/length(nmarker_values))
  teststat_vec[i] <- numerator/denominator
}
teststat_vec

#########

teststat_vec2 <- rep(NA, length(method_vec))
names(teststat_vec2) <- method_vec
for(i in 1:length(method_vec)){
  method_name <- method_vec[i]
  marker_values <- apply(glf_GPC[, names(glf_GPC)[grep(method_name, names(glf_GPC))] ], 2, median)
  nmarker_values <- apply(glf_nGPCblk[, names(glf_nGPCblk)[grep(method_name, names(glf_nGPCblk))] ], 2, median)
  
  numerator <- mean(marker_values) - mean(nmarker_values)
  denominator <- sqrt(var(marker_values)/length(marker_values) + var(nmarker_values)/length(nmarker_values))
  teststat_vec2[i] <- numerator/denominator
}
teststat_vec2

#########

teststat_vec3 <- rep(NA, length(method_vec))
names(teststat_vec3) <- method_vec
for(i in 1:length(method_vec)){
  method_name <- method_vec[i]
  marker_values <- apply(ery_Mark[, names(ery_Mark)[grep(method_name, names(ery_Mark))] ], 2, median)
  nmarker_values <- apply(ery_nMark[, names(ery_nMark)[grep(method_name, names(ery_nMark))] ], 2, median)
  
  numerator <- mean(marker_values) - mean(nmarker_values)
  denominator <- sqrt(var(marker_values)/length(marker_values) + var(nmarker_values)/length(nmarker_values))
  teststat_vec3[i] <- numerator/denominator
}
teststat_vec3


compute_t_test_vec <- function(data_Mark, data_nMark, f=median) {
  teststat_vec <- rep(NA, length(method_vec))
  names(teststat_vec) <- method_vec
  for(i in 1:length(method_vec)){
    method_name <- method_vec[i]
    marker_values <- apply(data_Mark[, names(data_Mark)[grep(method_name, names(data_Mark))] ], 2, f)
    nmarker_values <- apply(data_nMark[, names(data_nMark)[grep(method_name, names(data_nMark))] ], 2, f)
    
    numerator <- mean(marker_values) - mean(nmarker_values)
    denominator <- sqrt(var(marker_values)/length(marker_values) + var(nmarker_values)/length(nmarker_values))
    teststat_vec[i] <- numerator/denominator
  }
  names(teststat_vec)[4] <- 'velovi'
  teststat_vec
}

compute_t_test_vec(data_Mark=ery_Mark, data_nMark=ery_nMark, f=median)
compute_t_test_vec(data_Mark=ery_Mark, data_nMark=ery_nMark, f=mean) # mean of mean

compute_t_test_vec(data_Mark=glf_GPC, data_nMark=glf_nGPCblk, f=median)
compute_t_test_vec(data_Mark=glf_GPC, data_nMark=glf_nGPCblk, f=mean) # mean of mean

compute_t_test_vec(data_Mark=pan_Mark, data_nMark=pan_nMark, f=median) 
compute_t_test_vec(data_Mark=pan_Mark, data_nMark=pan_nMark, f=mean)  # mean of mean



apply(ery_Mark[, names(ery_Mark)[grep('utv', names(ery_Mark))] ], 2, median)
apply(ery_Mark[, names(ery_Mark)[grep('utv', names(ery_Mark))] ], 2, mean)


apply(ery_nMark[, names(ery_Mark)[grep('utv', names(ery_nMark))] ], 2, median)
apply(ery_nMark[, names(ery_Mark)[grep('utv', names(ery_nMark))] ], 2, mean)

par(mfrow=c(1,3))
compute_t_test_vec(data_Mark=pan_Mark, data_nMark=pan_nMark, f=mean) %>% barplot(main='pancreas',col='ivory')
compute_t_test_vec(data_Mark=ery_Mark, data_nMark=ery_nMark, f=mean) %>% barplot(main='erythroid',col='ivory')
compute_t_test_vec(data_Mark=glf_GPC, data_nMark=glf_nGPCblk, f=mean) %>% barplot(main='greenleaf',col='ivory')

compute_t_test_vec(data_Mark=glf_GPC, data_nMark=glf_nGPCblk, f=mean) %>% barplot(main='greenleaf',col='ivory')


# greenleaf
df <- compute_t_test_vec(data_Mark=glf_GPC, data_nMark=glf_nGPCblk, f=mean)
df <- data.frame( Method=names(df), Value=df )
df$Method <- factor(df$Method, levels = c("utv", "sct", "velovi_woprep", "velovi", "scv"))
# Barplot
(ggplot(df, aes(x = Method, y = Value)) + geom_bar(stat = "identity", fill = "lightsteelblue2") +
  labs(title = "", y = "Value", x = "Method") + theme_minimal() +
    theme( axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
           plot.title = element_text(size = 12, hjust = 0.5),  
           legend.text = element_text(size = 12), legend.title = element_text(size = 12) )) %>% 
  ggsave(filename='/Users/y2564li/Downloads/proj_scRNA/glf_barplot.png', width=6, height=4, bg='white')

# erythroid
df <- compute_t_test_vec(data_Mark=ery_Mark, data_nMark=ery_nMark, f=mean)
df <- data.frame( Method=names(df), Value=df )
df$Method <- factor(df$Method, levels = c("utv", "sct", "velovi_woprep", "velovi", "scv"))
# Barplot
(ggplot(df, aes(x = Method, y = Value)) + geom_bar(stat = "identity", fill = "lightsteelblue2") +
    labs(title = "Barplot of signal-to-random coherence score for erythroid dataset", y = "Value", x = "Method") + 
    theme_minimal() + theme( axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
                             plot.title = element_text(size = 12, hjust = 0.5),  
                             legend.text = element_text(size = 12), legend.title = element_text(size = 12) )) %>% 
  ggsave(filename='/Users/y2564li/Downloads/proj_scRNA/ery_barplot.png', width=6, height=4, bg='white')

# pancreas
df <- compute_t_test_vec(data_Mark=pan_Mark, data_nMark=pan_nMark, f=mean)
df <- data.frame( Method=names(df), Value=df )
df$Method <- factor(df$Method, levels = c("utv", "sct", "velovi_woprep", "velovi", "scv"))
# Barplot
(ggplot(df, aes(x = Method, y = Value)) + geom_bar(stat = "identity", fill = "lightsteelblue2") +
    labs(title = "Barplot of signal-to-random coherence score for pancreas dataset", y = "Value", x = "Method") + 
    theme_minimal() + theme( axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
                             plot.title = element_text(size = 12, hjust = 0.5),  
                             legend.text = element_text(size = 12), legend.title = element_text(size = 12) )) %>% 
  ggsave(filename='/Users/y2564li/Downloads/proj_scRNA/pan_barplot.png', width=6, height=4, bg='white')

