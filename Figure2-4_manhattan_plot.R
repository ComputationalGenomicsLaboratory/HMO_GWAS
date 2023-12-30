#This script was applied to generate the Manhattan plots presented in Main figures 1-4 and the supplementary figures of Ambalavanan et al.,

.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-pc-linux-gnu-library/3.4/") )
.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-pc-linux-gnu-library/3.5/") )
.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-redhat-linux-gnu-library/3.3/") )


#set your path in the HPC cluster 
setwd("/global/home/hpcxxxx/analysis/HMO/GWAS/all_HMOs")

#load your libraries
library(dplyr)
library(qqman)




#Plotting in loop
files <- list.files(path="/global/home/hpcxxxx/analysis/HMO/GWAS/all_HMOs", pattern=".assoc.linear", full.names = F, recursive = F)


for (i in 1:length(files)){
  res10 <- read.table(files[i], header = T, stringsAsFactors = F)
  res11 <- res10 %>% 
    filter(TEST == 'ADD')
  ylim <- abs(floor(log10(min(res11$P)))) + 2
  png(filename = paste(gsub("all_HMOs_multi\\.|INT\\_|.assoc.linear","", files[i]),'manhattan_multi','png', sep="."), width = 8, height = 4, units = 'in', res=300)
  manhattan(res11, logp=T, ylim=c(0, ylim),suggestiveline=-log10(1e-5), col = c("blue4","orange3"), cex=0.4, cex.axis = 0.6)
dev.off()
}

#Plot the figures in zoomed version
for (i in 1:length(files)){
  res10 <- read.table(files[i], header = T, stringsAsFactors = F)
  res11 <- res10 %>% 
    filter(TEST == 'ADD')
  #ylim <- abs(floor(log10(min(res11$P)))) + 2
  png(filename = paste(gsub("all_HMOs_multi\\.|INT\\_|.assoc.linear","", files[i]),'manhattan_multi_zoom','png', sep="."), width = 8, height = 4, units = 'in', res=300)
  manhattan(res11, logp=T, ylim=c(0, 15),suggestiveline=-log10(1e-5), col = c("blue4","orange3"), cex=0.4, cex.axis = 0.6)
dev.off()
}
####################################################END#######

