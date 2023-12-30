
#This script plots the overlaying manhattan plot of the meta analysis as shown in Figure 6

#library packages 
library(qqman)
library(dplyr)

#files of each HMO will be loaded in the loop
ldf<-list()
listfile <- dir(pattern = ".TBL")
for (k in 1:length(listfile)){
  ldf[[k]] <- read.table(listfile[k], header=T, stringsAsFactors = T) %>%
    mutate(HMO_Type=paste(gsub(".TBL","",listfile[k])))
}


#merging the 19 HMOs
all_hmo_19_suc <- do.call("rbind", ldf)

table(all_hmo_19_suc$HMO_Type)
#header <- scan("header", character(), quote = "")
#colnames(all_hmo_19_suc) <- c(header,"HMO_Type")

#extract the meta analysis results
meta.all.hmo <- which(all_hmo_19_suc[all_hmo_19_suc$Weight ==1375.00],) #total numbers including 395 and 980

meta.all.hmo <- all_hmo_19_suc[which(all_hmo_19_suc$Weight ==1375.00),]

write.table(meta.all.hmo, "Final_meta_analysis_results_hmo.txt", quote = F, row.names = F)


chr.data <- read.csv("../Sum.csv", header=T)
vars <- c("CHR", "SNP", "BP") 
subset.data <- chr.data[vars]

meta.data.chr <- merge(meta.all.hmo, subset.data, all.x = T, by.y = "SNP", by.x = "MarkerName")
write.table(meta.data.chr, "Final_meta_analysis_results_hmo2.txt", quote = F, row.names = F)

meta.data.chr <- cbind(meta.all.hmo, subset.data, by.)


######MANHATTAN PLOT##########
manhattan(meta.data.chr, snp="MarkerName" ,cex.axis = 0.8, col = c("#CC0621","#897074"), 
          main = "Meta analysis Manhattan Plot ", ylim=c(0,30), cex=0.8, p="P.value")

write.table(all_hmo_19_suc, "meta_analysis.txt", quote = F, row.names = F)
colnames(all_hmo_19_suc)

final_meta <- subset(meta.data.chr, select=c(1,9,10,6))












