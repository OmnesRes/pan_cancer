library(gplots)
intensities <- read.delim("average_coefficients.txt")
#-1 gives every column but the first
rmatrix1<-as.matrix(intensities[,-1])
row.names(rmatrix1) <- intensities[,1]


makeRowd <- function(intensities){
  rowcor <- as.dist(1-cor(t(intensities)))
  rowhc <- hclust(rowcor, method="average")
  rowd <- as.dendrogram(rowhc)
  return(list(rowd=rowd,rowhc=rowhc))
}


colors = c(seq(-3,-2,length=50),seq(-2,-1,length=50),seq(-1,1,length=100),seq(1,2,length=50),seq(2,3,length=50))
my_palette <- colorRampPalette(c("blue", "#4D4DFF","white", "#FF4D4D","red"))(n = 299)
pdf("clustergram.pdf", width=6,height=4)
tiff("clustergram_wkey.tiff",width=6,height=4,units="in",res=300)
heatmap.2(rmatrix1,
        Rowv=makeRowd(rmatrix1)$rowd,
        key=FALSE,
        #density.info=c("none"),
        #keysize=10,
        Colv=makeRowd(t(rmatrix1))$rowd,
        col=my_palette,
        breaks=colors,
        scale="col",
        trace="none",
        cexRow=.5,
        cexCol=1)

dev.off()


