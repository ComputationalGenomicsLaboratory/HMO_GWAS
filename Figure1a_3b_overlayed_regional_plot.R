#This script includes overlayed regional plot in Figure 2, 3 for chr 19 and 3

library(qqman)
library(ggplot2)
library(dplyr)

ldf<-list()
listfile <- dir(pattern = "assoc.linear")
for (k in 1:length(listfile)){
  ldf[[k]] <- read.table(listfile[k], header=F, stringsAsFactors = F) %>%
    mutate(HMO_Type=paste(gsub("top_sig\\.|all_HMOs_multi\\.|INT\\_|.assoc.linear","",listfile[k])))
}

all_hmo_19_suc <- do.call("rbind", ldf)
header <- scan("header", character(), quote = "")
colnames(all_hmo_19_suc) <- c(header,"HMO_Type")

manhattan(all_hmo_19_suc, cex.axis = 0.5, col = c("#CC0000", "#0000FF"), main = "Lactation GWAS Overlay Manhattan Plot (p<1e-5)", ylim=c(0,120))

CHR_19 <- all_hmo_19_suc %>% 
  filter(all_hmo_19_suc$CHR == 19) 

jpeg("regional_plot_chr19.jpg", width=6, height=4, units = "in", res=300)

ggplot(CHR_19, aes(x=CHR_19$BP, y=-log10(CHR_19$P), color= HMO_Type)) +
  geom_point(alpha = 1, size = 1.5) +
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



