# Corrinne Grover, November 2014
# Simple R script to calculate the percent identity values contained within an alignment
args <- commandArgs(TRUE)
FILE <- args[1]
OUT <- args[2]

library(seqinr)
library(ape)
alignment <- read.alignment(FILE, 'fasta', forceToLower = TRUE)
alignmentbin <- as.DNAbin(alignment)
bin <- dist.dna(alignmentbin, model = "JC", variance = FALSE, gamma = FALSE, pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE)
genemin <- 100*(1-min(bin))
genemax <- 100*(1-max(bin))
genemean <- 100*(1-mean(bin))
genemedian <- 100*(1-median(bin))

sink(OUT, append=TRUE)
cat(FILE)
cat("\t")
cat(genemin)
cat("\t")
cat(genemax)
cat("\t")
cat(genemean)
cat("\t")
cat(genemedian)
cat("\n")
sink()

