pkgs = c("knitr", "rmarkdown", "ggplot2", "ggpubr", "reshape2","cowplot",
         "superheat","plyr","dplyr","RColorBrewer")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("phyloseq")

