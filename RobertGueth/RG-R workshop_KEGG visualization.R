### Visualization of multiple RNAseq datasets in KEGG
### created by R Gueth on 11/02/15
### last updated: 11/02/15

# Overview - purpose:
# We'll take RNAseq data quantifying transcript abundances in treatment vs control for 1000's of genes
# and visualize them in various molecular pathways using KEGG
# Why? This is superior to the online visualization tool provided on the KEGG website since it 
# 1) enables much finer visualization control, 2) let's us visualize >1 dataset per pathway graph,
# and 3) allows for extensive and easy customization
#
# Overview - workflow:
# what we'll do is the following:
# 1. install all require packages and load the respective libraries
# 2. import some data and get it ready for visualization
# 3. retrieving KEGG pathway IDs
# 4. visualize the data in KEGG pathways

setwd("C:/Users/Robert/Desktop/Dropbox/stuff to work on/2015 Fall - BGSO R workshop")

###################################################
### 1. Packages & libraries ###
##############################
source("https://bioconductor.org/biocLite.R") # creates an instance of the newest Bioconductor installer used to d/l BioC packages
biocLite("KEGGREST") # installs the package 'KEGGREST' (case-sensitive) using the BioC installer biocLite
biocLite("pathview")
# packages need to be installed only once, but the corresponding libraries need to be initiated every time
library(KEGGREST)
library(pathview)


###################################################
### 2. Data import & prep ###
############################
# Import filtered trxome data table
trxome.data = read.table("RG_trxome filtered table 102714.csv",header=T,sep=",") # EntrezGene IDs in col 3
rownames(trxome.data) = trxome.data[,4]
trxome.data = trxome.data[,c(3,5:7)] # we'll only need the EntrezGene ID and expression value columns
trxome.data = as.matrix(trxome.data)   # for faster processing
# remove all genes with no EntrezGene ID using 'complete.cases'
gage.data = trxome.data[complete.cases(trxome.data[,1]),]
# remove all rows with duplicate EntrezGene IDs
dup.idx = which(duplicated(gage.data[,1]))
gage.data = gage.data[-dup.idx,]
# we'll change the rownames to EntrezGene IDs since 'pathview' needs it that way
rownames(gage.data) = gage.data[,1]
# we don't need the 1st column with EntrezGene IDs anymore, so let's remove it
gage.data = gage.data[,-1]


####################################################
### 3. Retrieving D. rerio KEGG pathway IDs ###
##############################################
### 3a. using KEGGREST ##
########################
kegg.list = keggList("pathway","dre") # keggList is a command from the library 'KEGGREST'
kegg.list = as.matrix(kegg.list)
path.ids = rownames(kegg.list)
# we want to retrieve the KEGG pathway IDs from the row names
# all row names have the format 'path:dre#####'
# we can do this in at least 3 different ways!
# 1: use regular expression to do this!
matches = regexec("path:dre([0-9]{5})",path.ids)
path.ids = regmatches(path.ids,matches)
path.ids = matrix(unlist(path.ids),ncol=2,byrow=T) # convert the list to a matrix with 2 columns row-wise for each list element
path.ids = path.ids[,2] # the pathway IDs are the items in the 2nd column
# 2: using the built-in string split function
path.ids = strsplit(x = path.ids,split = "dre",fixed = F) # returns a list with 2 items for each element
path.ids = matrix(unlist(path.ids),ncol=2,byrow=T) # convert the list to a matrix with 2 columns row-wise for each list element
path.ids = path.ids[,2] # the pathway IDs are the items in the 2nd column
# 3: using the improved string split function
biocLite("stringr") # the 'str_split_fixed' function is part of the 'stringr' package
library(stringr)
path.ids = str_split_fixed(string = path.ids,pattern = "dre",n = 2) # returns a matrix
path.ids = path.ids[,2] # the pathway IDs are the items in the 2nd column

######################
### 3b. using GAGE ##
####################
biocLite("gage")
library(gage)
kegg.dre = kegg.gsets(species="dre",id.type="kegg") # here we get all the KEGG pathways for D. rerio as a list
path.ids = attributes(kegg.dre$kg.sets) # this gets us the names of the elements on the list, which correspond to the pathway IDs and names
path.ids = matrix(unlist(path.ids),ncol=1,byrow=T) # convert the list to a matrix with 1 columns row-wise for each list element
# we'll use regular expression to get the pathway IDs
matches = regexec("dre([0-9]{5})",path.ids)
path.ids = regmatches(path.ids,matches) # we get a list with 2 items in each element, the 2nd of which we want
path.ids = matrix(unlist(path.ids),ncol=2,byrow=T) # convert the list to a matrix with 2 columns row-wise for each list element
path.ids = path.ids[,2] # the pathway IDs are the items in the 2nd column


###################################################
### 4. Visualization ###
#######################
# gage.data data columns:
# 1 - control tissue A, 2 - treatment1 tissue A, 3 - treatment2 tissue A
# setting some custom options for pathview:
## limit - fold change limit for up- and down-regulation (remember, this is in log2 scale)
## low, mid, high - color scale for down- to up-regulation
limit = list(gene=2,cpd=2) # this sets the upper and lower limit for the color scale on the pathway image (log2 scale)
low = list(gene="blue",cpd="blue") # this sets the color for the lower limit of the color scale
mid = list(gene="yellow",cpd="yellow") # this sets the color for the intermediate range of the color scale 
high = list(gene="red",cpd="red") # this sets the color for the upper limit of the color scale
gene.data.log2 = log2(as.numeric(gage.data[,2])/as.numeric(gage.data[,1])) # treatment1/control ratio, log2-transformed
gene.data.log2 = as.matrix(gene.data.log2)
rownames(gene.data.log2) = rownames(gage.data)
data = log2(as.numeric(gage.data[,3])/as.numeric(gage.data[,1])) # treatment2/control ratio, log2-transformed
gene.data.log2 = cbind(gene.data.log2,data) # appending the treatment2/control data to the data matrix
# In this last part we'll be creating the pathway images in the working directory
# we can either specify a single pathway or supply a vector of pathway IDs (here: path.ids)
# note that we are setting custom colors and color scale limits and are using log2-transformed data
# the node.sum option specifies how the expression value is determined for nodes that contain more than
# one gene (default is "sum", which tends to misrepresent what's actually happening)
# single pathway visualization
pv.out = pathview(gene.data.log2,pathway.id="04068",species="dre",limit=limit,low=low,mid=mid,high=high,node.sum="mean")
# all pathway visualization
# THIS WILL TAKE A WHILE!
pv.out = pathview(gene.data.log2,pathway.id=path.ids,species="dre",limit=limit,low=low,mid=mid,high=high,node.sum="mean")
