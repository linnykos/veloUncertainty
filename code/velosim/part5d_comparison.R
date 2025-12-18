rm(list=ls())

stochastic <- read.csv("~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/50-trials_stochastic.csv")[,1]
dynamical <- read.csv("~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/50-trials_dynamical.csv")[,1]

df <- cbind(stochastic, dynamical)
table(sign(df[,1] - df[,2]))

mean(df[,1]); sd(df[,1])
mean(df[,2]); sd(df[,2])