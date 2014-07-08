plot(itis[,1],xlim=c(0,10),type="b",col="blue",xlab="Coverage X10",ylab="No. of Insertion",main="No. of Insertion under different coverage",xaxp=c(0,10,10),yaxp=c(0,100,10),ylim=c(0,80))

lines(tif[,1],type="b",col="red")

legend("topleft",legend=c("ITIS","TIF"),col=c("blue","red"),lty=1)


