#This script includes creating the overlayed manhattan plot of conditional analysis and regional analysis plots 
#shown in Figure 4 (main figures)

library(qqman)
library(ggplot2)
library(dplyr)



#Loading the files of all the suggestive conditional analysis associations of plink output

setwd("~/SynologyDrive/CHILD/Projects_Queens/Lactation_all/HMO_paper_workspace_latest/Nature_Communications_revisions/Paper revisions/Analysis/Conditional_analysis/suggestive")

ldf<-list()
listfile <- dir(pattern = "assoc.linear")
for (k in 1:length(listfile)){
  ldf[[k]] <- read.table(listfile[k], header=F, stringsAsFactors = F) %>%
    mutate(HMO_Type=paste(gsub("suggestive\\.|conditional_analysis\\.|INT\\_|.assoc.linear","",listfile[k])))
}

all_hmo_19_suc <- do.call("rbind", ldf)
header <- scan("header", character(), quote = "")
colnames(all_hmo_19_suc) <- c(header,"HMO_Type")

write.table(all_hmo_19_suc, "../Conditional_analysis_all_hmo_suggestive_locus.txt", quote = F, row.names = F)

all_hmo_19_suc <- read.table("../Conditional_analysis_all_hmo_suggestive_locus.txt", header=T, sep="\t")
manhattan(all_hmo_19_suc, cex.axis = 0.8, col = c("#CC0621","#897074"), main = "Conditional GWAS Overlay Manhattan Plot (p<1e-5)", ylim=c(5,30), cex=0.8)


#Regional plot for chromosome 19
CHR_19 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 19) 

jpeg("../regional_plot_chr19.jpg", width=6, height=4, units = "in", res=300)

ggplot(CHR_19, aes(x=CHR_19$BP, y=-log10(CHR_19$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.0) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 19 (P-value<5e-5)") + 
  labs(x = "Chromosome (19p13.3 & 19q13.33)", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))


dev.off()

ggplot(CHR_19, aes(x=BP, y=BETA, color= HMO_Type)) +
  geom_point(alpha = 0.8, size = 1) +
  theme(legend.text=element_text(size=7)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0, size = 0.1, color = "red") +
  ggtitle("CHR 19") +
  labs(x = "Chromosome") +
  theme(plot.title = element_text(hjust = 0.5)) 

#Regional plot for chromosome 3

CHR_3 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 3) 

ggplot(CHR_3, aes(x=BP, y=BETA, color= HMO_Type)) +
  geom_point(alpha = 0.8, size = 1) +
  theme(legend.text=element_text(size=7)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0, size = 0.1, color = "red") +
  ggtitle("Chr 3") +
  labs(x = "Chromosome") +
  theme(plot.title = element_text(hjust = 0.5)) 

jpeg("../regional_plot_chr19.jpg", width=6, height=4, units = "in", res=300)


###############################Chr3 Regional Plot

jpeg("../regional_plot_chr3.jpg", width=6, height=4, units = "in", res=300)
ggplot(CHR_3, aes(x=CHR_3$BP, y=-log10(CHR_3$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 19 (P-value<1e-5)") + 
  labs(x = "Chromosome 3", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))

dev.off()


######################################Chr2 Regional Plot

CHR_2 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 2)

count.chr2 <- CHR_2 %>%
  filter(CHR_2$P<5e-08)

#unique=1

jpeg("../regional_plot_chr2.jpg", width=6, height=4, units = "in", res=300)
ggplot(CHR_2, aes(x=CHR_2$BP, y=-log10(CHR_2$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 2 (P-value<1e-5)") + 
  labs(x = "Chromosome 2", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))

dev.off()


#####################################Chr 7 Regional Plot

CHR_7 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 7)

count.chr7 <- CHR_7[which(CHR_7$P<5e-08),]
#unique=19

jpeg("../regional_plot_chr7.jpg", width=6, height=4, units = "in", res=300)
ggplot(CHR_7, aes(x=CHR_7$BP, y=-log10(CHR_7$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 7 (P-value<1e-5)") + 
  labs(x = "Chromosome 7", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))

dev.off()


########################################Chr16 regional Plot
CHR_16 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 16)

count.chr16 <- CHR_16[which(CHR_16$P<5e-08),]
#unique=7

jpeg("../regional_plot_chr16.jpg", width=6, height=4, units = "in", res=300)
ggplot(CHR_16, aes(x=CHR_16$BP, y=-log10(CHR_16$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 16 (P-value<1e-5)") + 
  labs(x = "Chromosome 16", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))

dev.off()
###########################################Chr18 Regional Plot
CHR_18 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 18)
count.chr18 <- CHR_18[which(CHR_18$P<5e-08),]
#unique=1



jpeg("../regional_plot_chr18.jpg", width=6, height=4, units = "in", res=300)
ggplot(CHR_18, aes(x=CHR_18$BP, y=-log10(CHR_18$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) +
  geom_hline(yintercept = -log10(5e-8), size = 0.1, color = "red") +
  ggtitle ("Chr 18 (P-value<1e-5)") + 
  labs(x = "Chromosome 18", y="-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size=10,colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size=10, colour = "black"))

dev.off()



