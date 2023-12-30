
#Load your library path

.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-pc-linux-gnu-library/3.4/") )
.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-pc-linux-gnu-library/3.5/") )
.libPaths( c( .libPaths(), "/global/home/hpc4094/R/x86_64-redhat-linux-gnu-library/3.3/") )


setwd("/global/home/hpc4094/analysis/HMO/GWAS/ceu_gwas2")

library(dplyr)
library(qqman)

#files<-list.files(path="/global/home/hpc4094/analysis/HMO/GWAS/all_HMOs", pattern=".assoc.linear", full.names = F, recursive = F)


##############################################################################################################
#Plotting in loop
files <- list.files(path="/global/home/hpc4094/analysis/HMO/GWAS/ceu_gwas2", pattern=".assoc.linear", full.names = F, recursive = F)


for (i in 1:length(files)){
  res10 <- read.table(files[i], header = T, stringsAsFactors = F)
  res11 <- res10 %>% 
    filter(TEST == 'ADD')
  ylim <- abs(floor(log10(min(res11$P)))) + 2
  png(filename = paste(gsub("all_HMOs_Ceu\\.|INT\\_|.assoc.linear","", files[i]),'manhattan_ceu','png', sep="."), width = 8, height = 4, units = 'in', res=600)
  manhattan(res11, logp=T, ylim=c(0, ylim),suggestiveline=-log10(1e-5), col = c("blue4", "cyan3"), cex=0.5, cex.axis = 1.0)
dev.off()
}

###################################################################################################################
#for (i in 1:length(files)){
 # res10 <- read.table(files[i], header = T, stringsAsFactors = F)
 # res11 <- res10 %>% 
 #   filter(TEST == 'ADD')
 # #ylim <- abs(floor(log10(min(res11$P)))) + 2
 # png(filename = paste(gsub("all_HMOs_multi\\.|INT\\_|.assoc.linear","", files[i]),'manhattan_multi_zoom','png', sep="."), width = 8, height = 4, units = 'in', res=300)
 # manhattan(res11, logp=T, ylim=c(0, 15),suggestiveline=-log10(1e-5), col = c("blue4", "orange3"), cex=0.4, cex.axis = 0.6)
#dev.off()
}
####################################################

#jpeg("DSLNH_multiethnic_QQ.jpeg") 

#qq(res11$P)
#dev.off()

#multiethnic_colors =c("#CC0621","#897074")
