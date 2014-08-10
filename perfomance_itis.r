svg("recoverate_on_simu.svg")


 itis <- read.table("itis_sort",as.is=T,row.names=1)
relocTE <- read.table("relocate_sort",as.is=T,row.names=1)
retroseq <- read.table("retro_sort",as.is=T,row.names=1)
temp <- read.table("TEMP_sort",as.is=T,row.names=1)
tif <- read.table("tif_sort",as.is=T,row.names=1)

total <- matrix(c(itis[,2],relocTE[,2],retroseq[,2],temp[,2],tif[,2]),ncol=5)

 overlap  <- matrix(c(itis[,3],relocTE[,3],retroseq[,3],temp[,3],tif[,3]),ncol=5)
exact_overlap  <- matrix(c(itis[,4],relocTE[,4],retroseq[,4],temp[,4],tif[,4]),ncol=5)

nm <- c("ITIS","RelocaTE","RetroSeq","TEMP","TIF")

bp <- barplot(total,beside=T,col="red",names.arg=nm,main="Recover Rate of Insertins Using Simulated Data")
barplot(overlap,beside=T,add=T,col="blue")
barplot(exact_overlap,beside=T,col="green",add=T)

 abline(h=52,lty="13")
text(2,54,labels="REF:52")
legend("topleft",legend=c("Total","Overlap 100bp","Exact overlap"),fill=c("red","blue","green"),bty="n")

mtext(side=1,cex=0.7,text= c("3x","5x","10x","15x","20x","25x"),at=bp,las=2)

dev.off()
