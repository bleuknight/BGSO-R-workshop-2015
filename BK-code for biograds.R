# load script that installs bioconductor
source("https://bioconductor.org/biocLite.R")
# Each of these commands tells Bioconductor to download and install each package

biocLite("org.Hs.eg.db" )
biocLite("DESeq2")
biocLite("parathyroidSE")
biocLite("GenomicRanges")
biocLite("GenomicFeatures")
biocLite("Rsamtools")
biocLite("DESeq2")
biocLite("lambda.r")
install.packages("Rcpp")
install.packages("acepack")


#access the sample BAM files
library("parathyroidSE")
extDataDir <- system.file("extdata", package = "parathyroidSE", mustWork = TRUE)
list.files( extDataDir )

#Experimental table containing metadata regarding conditions
sampleTable <- data.frame(
   sampleName = c( "Ctrl_24h_1", "Ctrl_48h_1", "DPN_24h_1" ),
   fileName = c( "SRR479052.bam", "SRR479053.bam",  "SRR479054.bam" ),
   treatment =  c( "Control", "Control", "DPN" ),
   time = c( "24h", "48h", "24h" ) )

#This is how the sample table should look
sampleTable

#construct the full paths to the bam files that we want to count
bamFiles <- file.path( extDataDir, sampleTable$fileName )
bamFiles

#We can peek into one of the BAM files to see the naming style of the 
#sequences (chromosomes). Here we use the BamFile function from the Rsamtools package.

library("Rsamtools")
seqinfo( BamFile( bamFiles[1] ) )

##We want to make sure that these sequence names are the same style as that 
##of the gene models we will obtain in the next section.

##You should have downloaded the GTF file before the start of class, make sure you
##put GTF file in your working directory (the same file where you save this code file)
getwd()
setwd ("C:/Users/Bleu/Documents/R")
library( "GenomicFeatures" )
hse <- makeTxDbFromGFF("C:/Users/Bleu/Documents/R/Homo_sapiens.GRCh37.75.subset.gtf", format="gtf" )

#Extract exons for each gene
exonsByGene <- exonsBy( hse, by="gene" )
exonsByGene

##We use the counting mode "Union", which indicates that those reads which 
##overlap any portion of exactly one feature are counted. As this experiment 
##produced paired-end reads, we specify singleEnd =FALSE. As protocol was not 
##strand-specific, we specify ignore.strand = TRUE. fragments = TRUE indicates 
##that we also want to count reads with unmapped pairs.
##This last argument is only for use with paired-end experiments.

library( "GenomicAlignments" )
se <- summarizeOverlaps( exonsByGene, BamFileList( bamFiles ), mode="Union",
                         singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )

##We investigate the resulting SummarizedExperiment class by looking at the 
##counts in the assay slot, the phenotypic data about the samples in colData 
##slot (in this case an empty DataFrame), and the data about the genes in the 
##rowData slot.

se

##check out count data, the "head" option just gives the first few lines of file

head( assay(se) )

#investigate column sum
colSums( assay(se) )

#Data about genes- all info about exons for each gene (ie each row)
rowRanges(se)

colData(se)

##assign meta data to colData

colData(se) <- DataFrame( sampleTable )

colData(se)$treatment

##We can extract columns from the colData using 
##the $ operator, and we can omit the colData to avoid extra keystrokes.
se$treatment


##We can also use the sampleName table to name the columns of our data matrix
colnames(se) <- sampleTable$sampleName
head( assay(se) )

str( metadata( rowRanges(se) ) )

##This summarized experiment is what we need to load into DeSeq2
# load the DESeq library
library(DESeq2)

#Run data
data( "parathyroidGenesSE" )
se <- parathyroidGenesSE
colnames(se) <- se$run

##Check column data
##Make sure that the order of rows in your column 
##data table matches the order of columns in the assay data slot.
colData(se)[1:5,1:4]



##A bonus about the workflow we have shown above is that information about
##the gene models we used is included without extra effort. The str 
##R function is used to compactly display the structure of the data in the list.
str( metadata( rowRanges(se) ) )


##Once we have our fully annotated SummerizedExperiment object, we can construct 
##a DESeqDataSet object from it, which will then form 
###the staring point of the actual DESeq2 package
ddsFull <- DESeqDataSet( se, design = ~ patient + treatment )

##Collapsing technical replicates is necessary because
##there are a number of samples which were sequenced in multiple runs.
##For example, sample SRS308873 was sequenced twice.
ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$sample,
                                    run = ddsFull$run )
head( as.data.frame( colData(ddsCollapsed)[ ,c("sample","runsCollapsed") ] ), 12 )

#as.data.frame forces R to show us the full list
head(as.data.frame( colData( ddsFull )[ ,c("sample","patient","treatment","time") ] ), 12)

##We can confirm that the counts for the new object are equal to the summed 
##up counts of the columns that had the same value for the grouping factor:
original <- rowSums( counts(ddsFull)[ , ddsFull$sample == "SRS308873" ] )
all( original == counts(ddsCollapsed)[ ,"SRS308873" ] )

##Here we will analyze a subset of the samples, 48 hours, with 
##either control, DPN or OHT treatment, 
##taking into account the multifactorial design.
dds <- ddsCollapsed[ , ddsCollapsed$time == "48h" ]

##It will be convenient to make sure that Control is the first level in the 
##treatment factor, so that the default log2 fold changes are calculated 
##as treatment over control and not the other way around. 
##The function relevel achieves this:
dds$treatment <- relevel( dds$treatment, "Control" )

##check and make sure samples are how we want them
as.data.frame( colData(dds) )

##remove genes that have expression levels of 0 for all samples
dds <- dds[ rowSums( counts(dds) ) > 0 , ]

##go back to specified design of experiment
design(dds)

#Run DESeq : Modeling counts with patient and treatment effects
dds <- DESeq(dds)
##A DESeqDataSet is returned which contains all the fitted 
##information within it, and the following section describes
##how to extract out results tables of interest from this object.



##Calling results without any arguments will extract the estimated log2 
##fold changes and p values for the last variable in the design formula. 
##If there are more than 2 levels for this variable (like now)
##results will extract the results table for a comparison of the 
##last level over the first level. 
res <- results( dds )
res

##view metadata for res
mcols(res, use.names=TRUE)

## view counts for gene with Smallest Pvalue
idx <- which.min(res$pvalue)
counts(dds)[idx, ]

##view normalized counts for same
counts(dds, normalized=TRUE)[ idx, ]

##get gene names
library( "org.Hs.eg.db" )

##the following function is for convenience and will circumvent the
##tedious process of using Annotationdbi- I'm not explaining it
#ids = list of IDS
#fromKey = key type; toKey = key type we want to convert to
#db = the AnnotationDb object to use.
#ifMultiple = the argument specifies what to do if one source ID maps to several target IDs:
#should the function return an NA or simply the first of the multiple IDs?
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
##END MYSTERY FUNCTION

res$hgnc_symbol <- convertIDs( row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db )
res$entrezid <- convertIDs( row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db )
head(res, 4)

##print out the results to file to work with in excel
write.table(res,file="DESeqOutput.txt",sep="\t")

##look at genes with adjusted P values below 0.1
sum( res$padj < 0.1, na.rm=TRUE )

##We subset the results table to these genes and then sort it by the log2 
##fold change estimate to get the 4 most down-regulated significant genes:
resSig <- subset(res, res$padj < 0.1 )
head( resSig[ order( resSig$log2FoldChange ), ], 4)

##and most upregulated
head( resSig[ order( -resSig$log2FoldChange ), ], 4)
