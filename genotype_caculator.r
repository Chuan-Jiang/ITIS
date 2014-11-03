#!/usr/bin/env  Rscript

opts <- as.numeric(commandArgs(trailingOnly = TRUE))
t <-binom.test(opts,p=0.5,alternative="greater")
#t <-binom.test(opts,p=2/3)
p <- t$p.value
cat(p)

