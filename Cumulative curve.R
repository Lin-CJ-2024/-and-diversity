
setwd("") 
rm(list = ls())
library(vegan)

otu <- read.delim('otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu=otu[,-dim(otu)[2]]
otu=t(otu)
sp <- specaccum(otu, method = 'random')
summary(sp)
plot(sp, ci.type = 'poly', col = 'blue', lwd = 0.1, ci.lty = 0.1, ci.col = 'black')
boxplot(sp, col = 'yellow', add = TRUE, pch = ' ')
