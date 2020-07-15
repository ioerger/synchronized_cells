#!/usr/bin/env Rscript

# usage: Rscript cluster.R <inputfile> <output_spreadsheet> <output_pdf> <K>
# input file: spreadsheet with curvefit means of expression for each each gene at each timepoint
# output file: tab-separated spreadsheet (.txt) with cluster numbers for each gene
# K: number of clusters (integer)
# example: Rscript cluster.R curvefit2_cos_means.txt curvefit2_cos_clust.txt curvefit2_cos_clusters.pdf 8

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outtxt = args[2]
outpdf = args[3]
K = args[4]

pdf(outpdf) # capture plots

data = read.table(infile,head=F,sep="\t")
rownames(data) = paste0(data[,1],"/",data[,2])
expr = data[4:18]
TP = c(3., 6.5, 9., 12., 18.5, 21., 27., 31., 33., 36., 39.5, 42., 45.5, 52., 55.)
colnames(expr) = paste0(TP,"hrs")
#plot(TP,expr[1,],type='l',xlab="Time (hrs)",ylab="Relative Expression")

d = dist(expr,method="euclidean") # euclidean is the default
d = dist(expr)
##hc = hclust(d,method="ward.D2")
hc = hclust(d,method="complete")

###############
#library(mclust)
#mod = mclustBIC(expr,G=1:25)
#plot(mod)
# I saved this as mclustBIC.png

# not very useful
# VVE peaks around 8; 
#Best BIC values:
#            VVE,8     VVE,7     VVE,9
#BIC      177.6242 102.31798  96.75646
###############

clust = cutree(hc,K)
table(clust)

temp = cbind(data[,1:3],clust,expr)
write.table(temp,outtxt,row.names=F,quote=F,sep='\t')

##par(mfrow=c(2,2))
par(mfrow=c(3,3))
for (k in 1:K)
{
  temp = expr[clust==k,]
  plot(TP,temp[1,],type='l',xlab="Time (hrs)",ylab="Relative Expression",main=paste("Cluster",k),ylim=c(-3,3))
  for (i in 2:length(temp[,1])) { 
    lines(TP,temp[i,],type='l',col=i) }
}

library(gplots)
heatmap.2(as.matrix(expr),trace="none",col=redblue(100),Colv=F,cexRow=0.3,keysize=1.2,lhei=c(1,4))

dev.off()
