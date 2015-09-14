library(gplots)



#-1 gives every column but the first
intensities <- read.delim("for_clustering.txt")
rmatrix1<-as.matrix(intensities[,-1])
row.names(rmatrix1) <- intensities[,1]




makeRowd <- function(intensities){
  rowcor <- as.dist(1-cor(t(intensities)))
  rowhc <- hclust(rowcor, method="average")
  rowd <- as.dendrogram(rowhc)
  return(list(rowd=rowd,rowhc=rowhc))
}


makeRowd(intensities)

colors = c(seq(-2.75,-2.5,length=25),seq(-2.5,-2,length=25),seq(-2,2,length=200),seq(2,2.5,length=25),seq(2.5,2.75,length=25))
my_palette <- colorRampPalette(c("blue", "#4D4DFF","white", "#FF4D4D","red"))(n = 299)
pdf("gene_clustering.pdf",, width=6,height=4)
tiff("gene_clustering.tiff",width=6,height=4,units="in",res=300)
heatmap.2(rmatrix1,
        Rowv=makeRowd(rmatrix1)$rowd,
        key=FALSE,
        #density.info=c("none"),
        #keysize=10,
        Colv=makeRowd(t(rmatrix1))$rowd,
        scale="row",
        col=my_palette,
        breaks=colors,
        trace="none",
        labRow=NA,
        cexCol=1,
        dendrogram="col")

dev.off()


