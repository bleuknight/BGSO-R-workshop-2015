#Read in a data table that we will use later

stress.data <- read.table("Data for Analysis.txt",head=T,sep="\t")

############## First, we are going to try this with randomly generated data that will show a pretty relationship #################3

# generates random data for acute cort measure
acute.cort <- rnorm(n=60, mean=102, sd=5.2)

#generates random data to be correlated with acute cort data
x <- acute.cort
correlatedValue = function(x, r){
  r2 = r**2
  ve = 1-r2
  SD = sqrt(ve)
  e  = rnorm(length(x), mean=0, sd=SD)
  y  = r*x + e
  return(y)
}

set.seed(5)
x = acute.cort
y = correlatedValue(x=x, r=.5)

#shows the correlation between the two variables
cor(x,y)

#re-labels the newly generated Y as the behavior score
behavior.score <- y


########################## Now we have random data, let's try the code! ##############################################

#We need to make a data frame with our species list to match our random data

#The species list comes from the stress.data that we read in at the beginning. Rename it to avoid confusion.
pretty.species <- stress.data$species

pretty.data <- data.frame(pretty.species, acute.cort, behavior.score)


# Installs the proper packages for phylogenetic analysis
install.packages("ape")
library("ape")

# Imports all 100 trees downloaded from birdtree.org
tree <- read.nexus("final.tree.nex",tree.names = NULL)
str(tree)

# Makes a consensus tree out of the 100 trees in the file "tree"
cons.tree <- consensus(tree,p=0.5,check.labels=TRUE)

# We wanted branch lengths for our tree, so this adds those
pretty.tree <- compute.brlen(cons.tree, method="Grafen", power=1)

#Take a look at the tree
plot(pretty.tree, cex=.6, font=1, adj=0, no.margin=TRUE)

#Need to make sure the tip labels match the species list in our data frame. 
pretty.tree$tip.label
pretty.data$pretty.species

#Uh oh, there is an underscore in the tip labels from the tree we downloaded, and just a space in our species list!

#Make a new vector called tree.tip.labels that consists of your tree tip labels
tree.tip.labels <- pretty.tree$tip.label

#And make a new vector out of the species names from your data frame
pretty.species.names <- pretty.data$pretty.species

#Use gsub to find a pattern, in this case the space between genus and species, and replace it with an underscore
pretty.species.names <- gsub(" ", "_", pretty.species.names)

#Replace the row names in your data frame with your newly formatted species names with the underscore
row.names(pretty.data)<-pretty.species.names


############################# Making the annotated phylogeny ##################################################################

# This installs the package necessary to run the annotated phylogeny code
install.packages("adephylo")
library(adephylo)

#Make a subset of only the data you want to look at
phylogeny.data <- subset(pretty.data, select = c(acute.cort, behavior.score))

# Makes an annotated phylogeny with trait data for acute cort and behavior score represented by a color heat map
pretty.annotated.phylogeny <- bullseye(pretty.tree, phylogeny.data, legend=FALSE, cex=.8, col.pal=list(function(x)rev(heat.colors(x))), 
                      traits.space=.04, traits.inset=1, label.offset=.1, rotate.tree=70, axis=FALSE)

# Can also use a different color for the second trait, but it doesn't always make patterns/correlations as clear
pretty.annotated.phylogeny.2 <- bullseye(pretty.tree, phylogeny.data, legend=FALSE, cex=.8, col.pal=list(function(x)rev(heat.colors(x)), function(y)rev(cm.colors(y))), 
                      traits.space=.04, traits.inset=1, label.offset=.1, rotate.tree=70, axis=FALSE)






