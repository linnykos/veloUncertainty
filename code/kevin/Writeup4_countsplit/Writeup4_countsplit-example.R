# https://anna-neufeld.github.io/countsplit.tutorials/articles/countsplit_tutorial.html
rm(list=ls())

library(countsplit)
library(Matrix)
library(ggplot2)
library(patchwork)

set.seed(1)
n <- 1000
p <- 200
X <- matrix(rpois(n*p, lambda=5), nrow=n)

set.seed(2)
split <- countsplit::countsplit(X, folds=2, epsilon=c(0.5,0.5))
Xtrain <- split[[1]]
Xtest <- split[[2]]

cor_vec <- sapply(1:p, function(j){
  stats::cor(Xtrain[,j], Xtest[,j])
})

hist(cor_vec)

cor_vec[1]
x <- Xtrain[,1] + stats::runif(n, min = -.25, max = .25)
y <- Xtest[,1] + stats::runif(n, min = -.25, max = .25)

png(file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup4/Writeup4_countsplit.png",
    height = 800, width = 800, units = "px", res = 300)
par(mar = c(4, 4, 0.5, 0.5))
plot(x, y, 
     asp = T, 
     col = rgb(0.5,0.5,0.5,0.1),
     pch = 16, 
     xlab = "Training",
     ylab = "Testing")
graphics.off()
