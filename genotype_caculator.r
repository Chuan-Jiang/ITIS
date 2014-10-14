#!/usr/bin/env  Rscript

opts <- as.numeric(commandArgs(trailingOnly = TRUE))
t <-binom.test(opts,p=2/3,alternative="greater")
p <- t$p.value
cat(p)

