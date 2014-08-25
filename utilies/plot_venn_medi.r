library("VennDiagram")
c <- readLines("count_ins_tools.tsv")
cl <- strsplit(c,split="\t")
nm <- lapply(cl,function(x){x <- x[1]})
cl <- lapply(cl,function(x){x <- x[-1]})
names(cl) <- nm

venn.diagram(cl,filename="venn.png",imagetype="png",fill=rainbow(5))
