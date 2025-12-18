rm(list=ls())

stochastic <- read.csv("~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/50-trials_stochastic.csv")[,1]
dynamical <- read.csv("~/kzlinlab/projects/veloUncertainty/out/kevin/simulation/50-trials_dynamical.csv")[,1]

df <- cbind(stochastic, dynamical)
table(sign(df[,1] - df[,2]))

mean(df[,1]); sd(df[,1])
mean(df[,2]); sd(df[,2])

#######

library(ggplot2)
library(tidyr)
library(ggpubr)

# 1. Convert your matrix/dataframe to a "long" format
# This creates a column for the model type and a column for the values
df_long <- as.data.frame(df) %>%
  pivot_longer(cols = everything(), 
               names_to = "Model", 
               values_to = "Value")

# 2. Create the boxplot with significance test
plot1 <- ggplot(df_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # Boxplot
  geom_jitter(width = 0.2, alpha = 0.5) +          # Optional: show individual points
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",          # Show actual p-value
                     label.x = 1.4) +              # Adjust position
  theme_minimal() +
  labs(title = "Comparison of Stochastic vs Dynamical Values",
       y = "Value",
       x = "Model Type") +
  scale_fill_brewer(palette = "Set2")

ggplot2::ggsave("~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/fig/kevin/scvelo/simulation_50-trials.png",
                plot1,
                height = 5, width = 7)
