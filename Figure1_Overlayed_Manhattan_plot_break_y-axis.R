#This is a submission script of the easy strata R package to obtain the overlayed manhattan plot of Figure 1, Ambalavanan et al.,
#Load your R library path and package
.libPaths( c( .libPaths(), "/global/home/hpcxxxx/R/x86_64-pc-linux-gnu-library/3.5/") )
library(EasyStrata)

#provide the full path of the ecf containing the script
EasyStrata("/global/home/hpcxxxx/analysis/HMO/GWAS/multi_GWAS/overlay_easystrata/MHplot_break_yaxis_HMO.ecf")

