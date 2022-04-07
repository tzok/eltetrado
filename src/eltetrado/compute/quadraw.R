#! /usr/bin/env Rscript

# source: http://colorbrewer2.org/?type=qualitative&scheme=Paired&n=9
o_plus      <- '#1f78b4'
o_minus     <- '#a6cee3'
n_plus      <- '#33a02c'
n_minus     <- '#b2df8a'
z_plus      <- '#ff7f00'
z_minus     <- '#fdbf6f'
default     <- 'gray'
canonical   <- 'black'

args = commandArgs(trailingOnly=TRUE)

# example data for testing
#
# seq.fasta:
# GUGUGUGU
#
# h1.helix:
# #8
# i j length value
# 1 3 1 1
# 5 7 1 1
#
# h2.helix
# #8
# 1 i length value
# 3 5 1 1
# 1 7 1 1
if (length(args) != 4) {
    stop("Usage: ./mp.R fasta helix1 helix2 pdf\nExample: ./mp.R /tmp/seq.fasta /tmp/h1.helix /tmp/h2.helix /tmp/output.pdf")
}


# install main dependencies
packages <- c('BiocManager')
check <- match(packages, installed.packages()[,1])
toinstall <- packages[is.na(check)]
if (length(toinstall) > 0) {
    install.packages(toinstall, repos = 'https://cloud.r-project.org/')
}

# install bioconductor dependencies
packages <- c('R4RNA')
check <- match(packages, installed.packages()[,1])
toinstall <- packages[is.na(check)]
if (length(toinstall) > 0) {
    BiocManager::install(toinstall, update = FALSE)
}

library(R4RNA)

fasta<-as.character(readBStringSet(args[1]))
helix1<-readHelix(args[2])
helix2<-readHelix(args[3])

helix1$col                          <- default
helix1$col[which(helix1$value==1)]  <- o_plus
helix1$col[which(helix1$value==2)]  <- o_minus
helix1$col[which(helix1$value==3)]  <- n_plus
helix1$col[which(helix1$value==4)]  <- n_minus
helix1$col[which(helix1$value==5)]  <- z_plus
helix1$col[which(helix1$value==6)]  <- z_minus
helix1$col[which(helix1$value==7)]  <- default
helix1$col[which(helix1$value==8)]  <- canonical

helix2$col                          <- default
helix2$col[which(helix2$value==1)]  <- o_plus
helix2$col[which(helix2$value==2)]  <- o_minus
helix2$col[which(helix2$value==3)]  <- n_plus
helix2$col[which(helix2$value==4)]  <- n_minus
helix2$col[which(helix2$value==5)]  <- z_plus
helix2$col[which(helix2$value==6)]  <- z_minus
helix1$col[which(helix2$value==7)]  <- default
helix1$col[which(helix2$value==8)]  <- canonical

x<-attributes(fasta)
x<-as.graphicsAnnot(x)
plotDoubleCovariance(helix1, helix2, top.msa=fasta[[1]], bot.msa=NA, add=FALSE, grid=TRUE, legend=FALSE, scale=FALSE, text=TRUE, lwd=3, pdf = args[4])

dim <- par("usr")
text(dim[1] + 1, dim[3] + 1, x$name[1], adj = 0, cex = 0.75, xpd = TRUE)
