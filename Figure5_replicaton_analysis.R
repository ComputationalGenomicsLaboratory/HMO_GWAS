
#This script includes the manhattan plot of replication analysis Figure 5
library(qqman)

idaho.final<- idaho %>% 
  filter(idaho$HMO == '2FL')
manhattan(idaho.final, chr="Chromosome_GRCh38", 
          bp="Position_GRCh38", 
          p="P.Value", 
          snp="rsIDs" )

manhattan(idaho, chr="Chromosome_GRCh38", 
          bp="Position_GRCh38", 
          p="P.Value", 
          snp="rsIDs" ,
          col=c("red", "pink1"), ylim=c(0,20), )


if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("ggmanh")

library(ggmanh)



manhattan_plot(x = idaho.final,
               pval.colname = "P.Value", 
               chr.colname = "Chromosome_GRCh38", 
               pos.colname = "Position_GRCh38", 
               plot.title = "Replication results", 
               rescale = TRUE, signif = c(5.6e-05, 2.2e-04), 
               highlight.col = highlight_colormap,
               color.by.highlight = TRUE)
idaho$color <-"Not Significant"
idaho$color[idaho$P.Value <= 2.2e-04] <- "Significant"

highlight_colormap <- c("Not Significant" = "grey", "Significant" = "blue")

tmp <- manhattan_data_preprocess(
  idaho, pval.colname = "P.Value", chr.colname = "Chromosome_GRCh38", pos.colname = "Position_GRCh38",
  highlight.colname = "color", highlight.col = highlight_colormap, signif = c(5.6e-05, 2.2e-04)
)
manhattan_plot(tmp, plot.title = "Replication Results",
               color.by.highlight = TRUE, rescale=T,
               rescale.ratio.threshold = 5)
