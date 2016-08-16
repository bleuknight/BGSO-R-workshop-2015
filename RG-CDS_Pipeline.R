### CDS extraction from contig sequences
### created by R Gueth on 11/02/15
### last updated: 11/02/15

# Overview - purpose:
# We've identified some contigs of interest in our transcriptome. Now we want to retrieve its
# coding sequence (CDS) in order to make sequence alignments, predict protein sequences, 
# construct phylogenetic trees, or design PCR primers;
# this pipeline will allow us to do so
#
# Overview - workflow:
# what we'll do is the following:
# 1. install necessary packages and run libraries
# 2. import contig sequences and local blastX results
# 3. execute the custom 'sequence.trimmer' function to trim the contig sequences to their CDSs
# 4. output the CDSs

setwd("C:/Users/Robert/Desktop/Dropbox/stuff to work on/2015 Fall - BGSO R workshop")

###################################################
### 1. Packages & libraries ###
##############################
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings") # the 'Biostrings' package includes functions to deal with nucleotide sequences
library(Biostrings)


###################################################
### 2. Data import & prep ###
############################
# importing file that was used for sequence submission to blast in order to get contig names and 
# sequences; make it into a 2-column matrix
cont.seq = read.table("Contig sequences for local blast.csv", header=F) # creates a list object
cont.seq = unlist(cont.seq)
# we want to order the contig names and corresponding sequence into the same row
# so we convert the character vector into a 2-column matrix
contig.ID.seq = matrix(cont.seq,ncol=2,byrow=T)
colnames(contig.ID.seq) = c("Contig ID","Contig sequence")

# import blastx results
blastx.results = read.table("Contig sequences from local blast.txt", header=F, sep="\t")
# we'll move columns 8 & 9 to the end of the list
blastx.results = cbind(blastx.results,blastx.results[,8:9])
blastx.results = blastx.results[-c(8:9)]
# let's put some column headers
# these are parameters that are returned by the local blastx search; we don't have to worry about
# what each of them represents except for the ones we need to use here:
# qseqid - query sequence ID, i.e. the contig name; qstart - query alignment start position (nucleotide);
# qend - query alignment end position (nucleotide); qframe - query reading frame for the alignment;
blastx.headers = c("qseqid", "sacc", "sallacc", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "pident", "qframe", "qseq", "sseq")
colnames(blastx.results) = blastx.headers
# The blastx search was set up to return the top 5 alignments for each contig sequence queried;
# We want to extract the information (qstart, qend, qframe) for the top hit only;
# so, we split the blastx data into sub-lists by unique contig names so that we can order
# the blast results for each contig query by bitscore
# to do that we use the 'split' function
blastx.subsets = split(blastx.results,blastx.results[,1]) # splits the blastx.results list into sub-lists by unique entries in the first column
attributes(blastx.subsets) # 13 elements/sub-lists, 1 for each input contig
# The following loop will be applied to each sublist in turn 
# For each sublist it will first sort the items in the sublist by decreasing value in the 'bitscore'
# Then it will extract the desired information: name (col1), frame (column 11 'qframe'), 
# start position (4 'qstart'), and end position (5 'qend')
info.collector = {} # we'll collect all information in here
contig.num = length(attributes(blastx.subsets)$names) # this gives us the number of contigs
for(i in seq(1,contig.num)){
  rank.bitscore = order(blastx.subsets[[i]]["bitscore"], decreasing=T) # order each sublist by bit score
  blast.info = blastx.subsets[[i]][rank.bitscore[1],c(1,4,5,11)] # extract the info we want for the first (=highest) bit score in each sublist
  # retrieve and add the contig sequence to info.collector as well
  contig.ID = paste(">",as.character(unlist(blast.info[1])),sep="")
  contig.seq = contig.ID.seq[match(contig.ID,table = contig.ID.seq[,1]),2] # looking up the contig sequence by contig.ID in the contig.ID.seq matrix
  blast.info = cbind(blast.info,contig.seq) # adding the contig sequence as an add'l column to blast.info
  info.collector = rbind(info.collector, blast.info) # for each contig we add blast.info as a new row to info.collector
}
# time to clean up the workspace!
rm(i,rank.bitscore,blast.info,contig.seq,contig.num,contig.ID,blastx.headers)
rownames(info.collector)=NULL


###################################################
### 3. Sequence trimming ###
###########################
# function sequence.trimmer: contig sequences are checked for orientation, extracted, trimmed,
# reverse-complemented, and collected in info.collector before being written as TXT output
# info.collector contains 5 variables: 1 - sequence ID, 2 - CDS start position, 3 - CDS stop position, 
# 4 - frame, 5 - contig sequence
sequence.trimmer = function(x){
  seq.trimmed = {} # we'll use 'seq.trimmed' to collect all trimmed sequences
  cat("Number of input contigs:",nrow(x)) # console output shows the number of sequences being processed - should be equal to the number of contigs we're processing
  for(i in seq(1,nrow(x))){
    if(x[i,4]<0){ # if the reading frame is negative
      extract.seq = toString(x[i,5]) # retrieve the contig sequence
      trim.seq = substring(extract.seq,x[i,3],x[i,2]) # use the start and end positions of the CDS alignment to trim the contig sequence; note how their order is flipped due to the way blastx returns them
      revcomp.seq = reverseComplement(DNAString(trim.seq))
      trim.seq = toString(revcomp.seq)
    }
    else(trim.seq=substring(toString(x[i,5]),x[i,2],x[i,3])) # for reading frames that are positive
    seq.trimmed=append(seq.trimmed,trim.seq)
  }
  return(seq.trimmed)
}

### Call sequence trimmer function
seq.trimmed = sequence.trimmer(info.collector)
# add the trimmed (CDS) sequences for each contig to the info.collector as a new column
info.collector = cbind(info.collector,seq.trimmed)
# get contig name and CDS sequence ready for export
trim.seq.out = unlist(t(info.collector[,c(1,6)])) # t() performs a transposition (rows<->columns)


###################################################
### 4. Output CDSs ###
#####################
write(trim.seq.out,"CDS-trimmed contig sequences.fa", sep="\n")
