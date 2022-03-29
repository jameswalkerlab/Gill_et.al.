library(tidyverse)
setwd('~/Dropbox/Broad-Research/ChemBio/R_Drpbx/sleep/alg10_opentargets_gwas_figure/')

v10<-read_delim('alg10_opentargets_variants.txt',delim='_',col_names=TRUE)
v10$xpos<-v10$pos-v10$pos[1]


v10b<-read_delim('alg10b_opentargets_variants.txt',delim='_',col_names=TRUE)
v10b$xpos<-v10b$pos-v10b$pos[1]
