#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Usage: ./mp.R fasta helix1 helix2 pdf\nExample: ./mp.R /tmp/seq.fasta /tmp/h1.helix /tmp/h2.helix /tmp/output.pdf")
}

library(R4RNA)
library(Biostrings)
library(RColorBrewer)

fasta<-as.character(readBStringSet(args[1]))
helix1<-readHelix(args[2])
helix2<-readHelix(args[3])

helix1$col<- "green"
helix1$col[which(helix1$value==1)] <- "black"
helix1$col[which(helix1$value==2)] <- "gray50"
helix1$col[which(helix1$value==3)] <- "blue"
helix1$col[which(helix1$value==4)] <- "cyan"
helix1$col[which(helix1$value==5)] <- "red"
helix1$col[which(helix1$value==6)] <- "orange"

helix2$col<- "green"
helix2$col[which(helix2$value==1)] <- "black"
helix2$col[which(helix2$value==2)] <- "gray50"
helix2$col[which(helix2$value==3)] <- "blue"
helix2$col[which(helix2$value==4)] <- "cyan"
helix2$col[which(helix2$value==5)] <- "red"
helix2$col[which(helix2$value==6)] <- "orange"

x<-attributes(fasta)
x<-as.graphicsAnnot(x)
plotDoubleCovariance(helix1, helix2, top.msa=fasta[[1]], bot.msa=NA, add=FALSE, grid=TRUE, legend=FALSE, scale=FALSE, text=TRUE, lwd=3, pdf = args[4])

dim <- par("usr")
text(dim[1] + 1, dim[3] + 1, x$name[1], adj = 0, cex = 0.75, xpd = TRUE)
