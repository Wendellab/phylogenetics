# an Rscript to calculate the number of variable positions over an entire alignment
# Corrinne Grover, December 2014

args <- commandArgs(TRUE)
FILE <- args[1]
OUT <- args[2]

library(seqinr)
library(ape)
alignment <- read.fasta(file = FILE, seqtype = "DNA", forceDNAtolower = TRUE)
Seq1 <- alignment[[1]]
Seq2 <- alignment[[2]]
Seq3 <- alignment[[3]]
Seq4 <- alignment[[4]]
Seq5 <- alignment[[5]]
Seq6 <- alignment[[6]]
Seq7 <- alignment[[7]]

varsites <- which(Seq1!=Seq2)
varsites <- append(varsites, which(Seq1!=Seq3))
varsites <- append(varsites, which(Seq1!=Seq4))
varsites <- append(varsites, which(Seq1!=Seq5))
varsites <- append(varsites, which(Seq1!=Seq6))
varsites <- append(varsites, which(Seq1!=Seq7))
varsites <- append(varsites, which(Seq2!=Seq3))
varsites <- append(varsites, which(Seq2!=Seq4))
varsites <- append(varsites, which(Seq2!=Seq5))
varsites <- append(varsites, which(Seq2!=Seq6))
varsites <- append(varsites, which(Seq2!=Seq7))
varsites <- append(varsites, which(Seq3!=Seq4))
varsites <- append(varsites, which(Seq3!=Seq5))
varsites <- append(varsites, which(Seq3!=Seq6))
varsites <- append(varsites, which(Seq3!=Seq7))
varsites <- append(varsites, which(Seq4!=Seq5))
varsites <- append(varsites, which(Seq4!=Seq6))
varsites <- append(varsites, which(Seq4!=Seq7))
varsites <- append(varsites, which(Seq5!=Seq6))
varsites <- append(varsites, which(Seq5!=Seq7))
varsites <- append(varsites, which(Seq6!=Seq7))

numvarsites = length(unique(varsites))
totalsites = length(Seq1)

percvarsites = 100*(numvarsites/totalsites)

sink(OUT, append=TRUE)
cat(FILE)
cat("\t")
cat(percvarsites)
cat("\n")
sink()
