library("VennDiagram")
c <- readLines("count_ins_tools.tsv")
cl <- strsplit(c,split="\t")
nm <- lapply(cl,function(x){x <- x[1]})
cl <- lapply(cl,function(x){x <- x[-1]})
names(cl) <- nm



#venn.diagram(cl,cat.cex=0.6,margin=0.2,filename="venn.png",imagetype="png",fill=rainbow(5))
venn.diagram(cl,cat.cex=0.7,cat.dist=c(0.1,0.08,0.1,0.1,0.08),alpha = 0.50,margin=0.2,filename="venn.png",imagetype="png",fill=rainbow(5), cat.fontface = "bold")
