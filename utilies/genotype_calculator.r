paras <- as.numeric(commandArgs(trailingOnly =  T))

p <- paras[3]
num <- paras[c(1,2)]

test <- binom.test(num,p=p)
pv <- test$p.value

Q <- -10*log10(pv)

cat(Q)
