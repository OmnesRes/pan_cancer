library(VennDiagram)

LIHC<-read.table('LIHC_bad.txt')
LIHC<-as.character(unlist(LIHC))
LUAD<-read.table('LUAD_bad.txt')
LUAD<-as.character(unlist(LUAD))
KIRP<-read.table('KIRP_bad.txt')
KIRP<-as.character(unlist(KIRP))

COAD<-read.table('COAD_good.txt')
COAD<-as.character(unlist(COAD))
GBM<-read.table('GBM_good.txt')
GBM<-as.character(unlist(GBM))
LUSC<-read.table('LUSC_good.txt')
LUSC<-as.character(unlist(LUSC))


venn.plot<-venn.diagram(list(A=LIHC, B=LUAD, C=KIRP),category.names=c('LIHC', 'LUAD', 'KIRP'),fill = c("#B23212", "#00FF00",'#00FFFF'),alpha = c(.9, .5,.5), cex = 2,cat.fontface = 1,cat.cex=c(2,2,2),lty =2,cat.dist=c(.06,.06,.04),filename=NULL)
                        #filename='lihc_luad_kirp_bad.tiff')

pdf('venn1.pdf')
grid.draw(venn.plot)
dev.off()

venn.plot<-venn.diagram(list(A=COAD, B=GBM, C=LUSC),category.names=c('COAD', 'GBM', 'LUSC'),fill = c("#B23212", "#00FF00",'#00FFFF'),alpha = c(.9, .5,.5), cex = 2,cat.fontface = 1,cat.cex=c(2,2,2),lty =2,cat.dist=c(.06,.06,.04),filename=NULL)



pdf('venn2.pdf')
grid.draw(venn.plot)
dev.off()
