# Heatmap and clustering based on log transformed read counts
# ============================================================
args = commandArgs(trailingOnly=TRUE)
library(gplots)
setwd(args[1])
in_filename = args[2]
delta = as.double(args[3])
out_filename = args[4]
sgcounts_File = read.csv(file=in_filename,head=TRUE,sep="\t")
sgcounts = sgcounts_File[1:nrow(sgcounts_File),3:(ncol(sgcounts_File)-1)]
log_sgcounts = log(sgcounts+delta)
DataMatrix = as.matrix(log_sgcounts)
col_palette = colorRampPalette(c("red","yellow","green"))(n = 40)
pdf(out_filename)
heatmap.2(DataMatrix,trace="none",col=col_palette,dendrogram="column",margins=c(7,7),labRow="")
dev.off()