# Heatmap and clustering based on log transformed read counts
# ============================================================

# Load heatmap.2
library(gplots)

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
in_filename = args[2]
delta = as.double(args[3])
out_filename = args[4]
W = as.integer(args[5])
H = as.integer(args[6])
S = as.integer(args[7])
r = as.integer(args[8])
svg_out = args[9]

# Read data
sgcounts_File = read.csv(file=in_filename,head=TRUE,sep="\t")
sgcounts = sgcounts_File[1:nrow(sgcounts_File),3:(ncol(sgcounts_File)-1)]

# Process data and plot
log_sgcounts = log(sgcounts+delta)
DataMatrix = as.matrix(log_sgcounts)
col_palette = colorRampPalette(c("yellow","orange","red"))(n = 40)
png(paste(out_filename,".png",sep=""), width = W, height = H, pointsize = S)
heatmap.2(DataMatrix,trace="none",col=col_palette,dendrogram="column",margins=c(r,r),labRow="")
dev.off()
if (svg_out == "True") {
  svg(paste(out_filename,".svg",sep=""))
  heatmap.2(DataMatrix,trace="none",col=col_palette,dendrogram="column",margins=c(r,r),labRow="")
  dev.off()
}