itis <- read.table("nf54_itis_counts.txt",as.is=T,row.names=2)
tif <- read.table("nf54_tif_counts.txt",as.is=T,row.names=2)
retro <- read.table("nf54_retro_counts.txt",as.is=T,row.names=2)
reloc <- read.table("nf54_reloca_counts.txt",as.is=T,row.names=2)


xpos = c(3,5,10,15,20,25,30,40,50,60,70,80,90,100)

plot(x=xpos,y=itis[,1],xlim=c(0,100),type="b",col="blue",xlab="Coverage",ylab="No. of Insertion",main="No. of Insertion under different coverage",xaxp=c(0,100,10),yaxp=c(0,100,10),axes=F)

lines(x=xpos,y=tif[,1],type="b",col="red")
lines(x=xpos,y=retro[,1],type="b",col="green")
lines(x=xpos,y=reloc[,1],type="b",col="black")
box()
axis(side=1,at=xpos,labels=F,tcl=-0.3)
axis(side=2,at=seq(10,100,by=10))

text(xpos,par("usr")[3]-1,srt=-45,adj=0,labels=paste(xpos,rep("x",14)),xpd=T,cex=0.8)

legend("topleft",legend=c("ITIS","TIF","retroSeq","relocaTE"),col=c("blue","red","'green","black"),lty=1)


