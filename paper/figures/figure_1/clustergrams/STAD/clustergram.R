##load the necessary libraries
library(gplots)
library(RColorBrewer)


##read the data
intensities <- read.delim("for_clustering_100_good_100_bad.txt")
rmatrix2<-as.matrix(intensities[,-1])
row.names(rmatrix2) <- intensities[,1]


##a distance function that uses pearson correlations
makeRowd <- function(intensities){
  rowcor <- as.dist(1-cor(t(intensities)))
  rowhc <- hclust(rowcor, method="average")
  print(cor(cophenetic(rowhc), rowcor))
  rowd <- as.dendrogram(rowhc)
  return(list(rowd=rowd,rowhc=rowhc))
}


makeRowd(intensities)

##want a sidebar
cl = read.table("genes_with_prognosis.txt", sep='\t')

color.map <- function(cl) {
  if(cl == "Good"){ 
    cl <- "green"
  }else if(cl == "Bad"){ 
    cl <- "red"
  }
  return(cl)
}

sidebarcolors <- apply(cl[2,],2, color.map)


##the tiff is for the final figure, the pdf is used when the patients in the clusters need to be identified
tiff("clustergram.tiff",width=6,height=4,units="in",res=300)
#pdf("clustergram_with_patients.pdf",width=10,height=4)


heatmap.2(rmatrix2,
          Rowv=FALSE,
          Colv=makeRowd(t(rmatrix2))$rowd,
          key=FALSE,
          #keysize=10,
          density.info=c("none"),
          col=colorRampPalette(c("yellow","blue"))(64),
          scale="row",
          trace="none",
          cexCol=.1,
          labRow=NA,
          ##when making the pdf comment out the labCol=NA
          labCol=NA,
          RowSideColors=sidebarcolors,
          dendrogram="none")

dev.off()

##writes the new order of the patients
write(makeRowd(t(rmatrix2))$rowhc$order,'column_order.txt')
