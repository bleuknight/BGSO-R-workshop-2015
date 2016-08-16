# Example 1
# Pattern description: genename@chromosome. gene name is in the format of 
#   one upper case letter+ lower case letters + number, 
#   chromosome is in the format of chrX, X can be 1-23, X, Y
gene1 <- "Rab21@chr1"
gene2 <- "Domo32@chr23"
gene3 <- "Assm2@chrY"

# pattern <- "(.*)@(.*)"
pattern <- "([A-Z][a-z]*[0-9]*)@chr[0-9]*[X,Y]*"
matched <- regexec(pattern, gene1)
values <- regmatches(gene1, matched)

gene <- c("Rab21@chr1", "Domo32@chr23", "Assm2@chrY")
pattern <- "([A-Z][a-z]*[0-9]*)@chr[0-9]*[X,Y]*"
matched <- regexec(pattern, gene)
values <- regmatches(gene, matched)
names <- c(values[[1]][2],values[[2]][2],values[[3]][2])



# Example mosquito genome
gff3 <- read.table(file = "Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.3.gff3.gz", nrows = 1000, sep = '\t')
# Objective: Find all gene names within gff3 annotation result.
#select all genes in my genome (in column 3)
tf.gene <- gff3$V3 == "gene"
#return all gene names
names <- gff3$V9[tf.gene]


# using regular expression to extract names.
pattern <- "ID=([^;]*)"
matched <- regexec(pattern, as.character(names))
values <- regmatches(as.character(names), matched)
all.names <- sapply(values, function(x) { x[2] })