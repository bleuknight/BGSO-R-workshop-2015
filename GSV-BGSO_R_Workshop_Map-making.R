# G.S.Vidaurre 3-Nov-15

###### Downloading and mapping GBIF data

# to download data from GBIF, you need to create an account online
# for the purposes of the workshop, I have commented the code to download GBIF data
# and provided you with the downloaded data, which you can read into RStudio

# install the Global Biodiversity Information Facility (GBIF) package
# install.packages("rgbif") # comment this line once package installed
# library(rgbif)

# start by obtaining taxon key for species of interest
# key <- name_suggest(q='Myiopsitta monachus', rank='species')$key[1]
# key # 2479407

# how many Myiopsitta monachus obsvervations have coordinates?
# occ_search(taxonKey=key, limit=20, hasCoordinate = TRUE) # 32636!

##### Downloading data from GBIF ######
# mopa <- occ_download('taxonKey = 2479407', user = "gsmithvi", pwd = "moParak33t",
#                      email = "gsmithvi@gmail.com")

# # check download status
# occ_download_meta(mopa)

# # when done running, download data
# mopa <- occ_download_get(mopa, overwrite = TRUE)

# # import into R and start manipulating
# # some warning messages about changing class of some data, mostly time of day
# suppressPackageStartupMessages(library("dplyr"))
# df <- tbl_df(occ_download_import(mopa))

# # saving data to Desktop to avoid having to download again
# setwd("~/Desktop")
# write.csv(df, "Myiopsitta_monachus_GBIF_3nov15.csv")
########################################

##### Reading GBIF data into RStudio #####
setwd("~/Desktop/BGSO_R_workshop/Map-making_Grace/")
mopa <- read.csv("Myiopsitta_monachus_GBIF_3nov15.csv", header = TRUE)

# checking data
names(mopa) # column namws
str(mopa) # information about data type per column
# View(mopa) # viewing data frame in new RStudio tab
unique(mopa$vernacularName) # are all these observations monk parakeets?
unique(mopa$scientificName) # the grebe doesn't show up by Latin name...

# filter out the single grebe observation
table(mopa$vernacularName)
mopa1 <- mopa[grep("grebe", mopa$vernacularName, ignore.case = TRUE, invert = TRUE), ]
nrow(mopa) 
nrow(mopa1) # new data frame has 1 less row

# 3238 observations have no common name
# do these have the correct scientific name?
mopa2 <- mopa1[grep("monk parakeet", mopa1$vernacularName, ignore.case = TRUE, invert = TRUE), 
               grep("scientificName$", names(mopa1))]
# we will keep these observations, as they all pertain to the Myiopsitta genus
length(mopa2[grepl("Myiopsitta", mopa2, ignore.case = TRUE)]) == length(mopa2)

# remove allospecies Myiopsitta luschi
length(grep("luchsi", mopa1$scientificName)) # 41 rows should be removed
mopa3 <- mopa1[grep("luchsi", mopa1$scientificName, invert = TRUE), ]
nrow(mopa1) - nrow(mopa3)

# now we will manipulate the data to retain only the columns we want
mopa4 <- mopa3[, grep(paste(c("decimalLatitude", "decimalLongitude", "^genus$", "speciesKey", "vernacularName" , "year", 
                         "taxonID", "scientificName$", "verbatimLocality", "stateProvince", 
                         "gbifID"), 
                       collapse = "|"), names(mopa3))]
str(mopa4)
names(mopa4)
nrow(mopa4)

# unfortunately, not much data available on which organizations contributed data
# grep("source", names(mopa1))
# names(mopa1)[50]
# table(mopa1$ownerInstitutionCode)
# table(mopa1$rights)
# table(mopa1$publisher)
# table(mopa1$license)

# remove observations without coordinates (e.g., with NAs)
length(which(is.na(mopa4$decimalLatitude) & is.na(mopa4$decimalLongitude))) # 775 should be removed
mopa5 <- mopa4[which(!is.na(mopa3$decimalLatitude) & !is.na(mopa4$decimalLongitude)), ]
nrow(mopa4) - nrow(mopa5)


##### Using maps package #####

# install.packages("maps")
library(maps)

lat <- mopa5$decimalLatitude
lon <- mopa5$decimalLongitude
ext <- 5 

# reinitialize after creating mopa6
# lat <- mopa6$decimalLatitude
# lon <- mopa6$decimalLongitude
# lat2 <- mopa5$decimalLatitude # place Antarctica back onto map

# divide lat and lon by invasive / native status
lat.inv <- lat[lat >= 0]
lon.inv <- lon[which(lat >= 0)]
lat.nat <- lat[lat <= 0]
lon.nat <- lon[which(lat <= 0)]

# creating an image file can be mroe useful than exporting plots from RStudio
# can change resolution for presentations/pubs
tiff(file = "Myiopsitta_monachus_GBIFmap.tiff", units = "cm", height = 10,
     width = 10, res = 300)

# create an empty map after the first time running code below,
# as well as every time you run dev.off()
map("world", xlim = c(min(lon) - ext, max(lon) + ext),
    ylim = c(min(lat2) - ext, max(lat2) + ext), interior = FALSE, fill = FALSE)

# par("usr") specifies plot boundaries using the limits of the data itself
# here we create a gridded rectangle as the background for the map
rect(par("usr")[1], par("usr")[3], par("usr")[2], 
         par("usr")[4], col = cm.colors(20)[7])
    abline(h = seq(-90, 90, 5), col = "white", lwd = 0.9)
    abline(h = seq(-90, 90, 5), col = "white", lwd = 1.1)
    abline(v = seq(-180, 180, 5), col = "white", 
           lwd = 0.9)
    abline(v = seq(-180, 180, 5), col = "white", 
           lwd = 1.1)
    
map("world", xlim = c(min(lon) - ext, max(lon) + ext), add = TRUE,
    ylim = c(min(lat2) - ext, max(lat2) + ext), fill = TRUE,
    col = terrain.colors(10)[6]) 

mtext("Longitude (DD)", side = 1, line = 1)
mtext("Latitude (DD)", side = 2, line = 1)
title("Myiopsitta monachus sightings", cex = 0.7, line = 1)

# using alpha argument to change transparency for better density visualization
# points(lon, lat, pch = 21, cex = 1, col = heat.colors(10)[2], 
#        bg = heat.colors(20, alpha = 0.2)[3])

# color points by native or invasive status
points(lon.inv, lat.inv, pch = 21, cex = 1, col = "black", lwd = 0.5,
       bg = heat.colors(20, alpha = 0.2)[4])

points(lon.nat, lat.nat, pch = 21, cex = 1, col = "black", lwd = 0.5,
       bg = topo.colors(10, alpha = 0.2)[3])

# add line for equator
abline(h = 0, col = "black", lwd = 1.1)

# add scale to map
map.scale(ratio = FALSE, relwidth = 0.2, x = 0, y = -40, cex = 0.5)

# turn off graphics device if seeing weird plot errors in RStudio
# always run dev.off() if creating tiff file
dev.off() 

# what are those funny points down by Antarctica?
# the stateProvince shows up as Salta or Chaco
# these must be data entry errors, so will remove them
# then rerun the mapmaking code with mopa6
# View(mopa5[which(mopa5$decimalLatitude < -50),])
mopa6 <- mopa5[!mopa5$decimalLatitude < -50,]
length(which(mopa5$decimalLatitude < -50)) == nrow(mopa5) - nrow(mopa6)


##### Using ggmaps package #####

# install.packages("ggmap")
library(ggmap)

# we will start with invasive populations in Florida
table(mopa6$stateProvince)
mopa6_FL <- mopa6[mopa6$stateProvince == "Florida", ]
nrow(mopa6_FL)
# View(mopa6_FL)

# remove all observations with latitude < 25 degrees
# from playing around, seems other countries have a province with the same name
mopa6_FL <- mopa6_FL[mopa6_FL$decimalLatitude >= 25, ]
nrow(mopa6_FL)

# I want to plot international US Customs ports in FL
# http://www.portcodes.com/ports.php?name=&state=FL&portcode=
# copied table (now table on clipboard)
# the code below imports the table on my clipboard into R (for Macs)
# FL.ports <- read.table(pipe("pbpaste"), sep = "\t", col.names = c("Number","Port", "State", "PortCode"),
#                        row.names = 1)

# pipe("xclip -i", "w") # UNIX 
# file(description = "clipboard"...) # Windows

# write.csv(FL.ports, "FL_ports.csv")
FL.ports <- read.csv("FL_ports.csv", header = TRUE, row.names = 1)
nrow(FL.ports)
str(FL.ports)

# now I can use the geocode function to download coordinates of these ports,
# as well as the capital of FL
gc <- geocode(location = c(as.character(FL.ports$Port), "Miami"), source = "google",
              messaging = FALSE)

# combine these geocode coords with Mymon data
# ggplot only allows aes to be called once, so we need all points in a single dataset

# create new Mymon Florida dataset with only lat and lon values
mopa7_FL <- data.frame(lon = mopa6_FL$decimalLongitude, lat = mopa6_FL$decimalLatitude)

# adding a "group variable" that will have levels for parakeet sightings, ports
# and the city of Miami, to use for ggplotting colors 
mopa7_FL_gc <- cbind(rbind(mopa7_FL, gc), group = c(rep(1, nrow(mopa6_FL)), 
                                                    rep(2, nrow(gc)-1), 3))
nrow(mopa7_FL_gc) - nrow(mopa7_FL) # should be 26

# convert group to factpr variable for color plotting purposes
mopa7_FL_gc$group <- as.factor(mopa7_FL_gc$group)

# change names of levels to have an appropriate legend 
levels(mopa7_FL_gc$group) <- c("mopa_sighting", "US Customs port", "Miami")
str(mopa7_FL_gc)

# download the map raster 
google <- get_googlemap("florida", zoom = 7)

# when creating the map, you can specify whether or not you want to change marker colors,
# sizes, fills by group (any factor variable in your dataset, here it is merely "group")
p <- ggmap(google) + geom_point(aes(x = lon, y = lat, color = group, fill = group, 
                                    size = group, shape = group), data = mopa7_FL_gc)  

# when adding layers, you can specify the colors, marker types, etc. for each factor level
# I'm initializing colors and fill separately to make the code less cluttered
cols <- c(heat.colors(20, alpha = 0.2)[4], "slateblue1", heat.colors(10)[9])
fill.cols <- rep(c(heat.colors(20, alpha = 0.2)[4], "slateblue1", "black"))

p + scale_fill_manual(values = cols) + scale_colour_manual(values= fill.cols) + 
  scale_shape_manual(values=c(21,24,22)) + scale_size_manual(values=c(2,3,4))
 

##### Using shapefiles with ggmaps #####

# now I want to add some information about FL human population density
# searched "urban" on http://www.fgdl.org/metadataexplorer/explorer.jsp
# downloaded "2010 U.S. CENSUS URBAN AREAS AND CLUSTERS IN FLORIDA"
# shows up as folder called "ua2010", which I moved to my Desktop

# NOTE: the code below reqires the package rgdal
# if this doesn't install on your first try using install.packages("rgdal"),
# then please hold back on running the code below, and instead watch as I run it
# I'm happy to help install the relevant packages on your computer if necessary

# install.packages("rgdal")
# install.packages("maptools")
library(sp)
library(rgdal)
library(maptools)

setwd("~/Desktop/BGSO_R_workshop/Map-making_Grace/ua2010/")
FLurb <- readOGR(".", "ua2010") # reading files for urban clusters

# transforming between datum and projection(s) using rgal in conjunction with PROJ.4
# e.g., transforming lcoation of points from curved surface to flat surface
FLurb <- spTransform(FLurb, CRS("+proj=longlat +datum=WGS84"))
# can also try subsituting this CRS("+proj=longlat +no_defs")

# now convert to a dataframe that R can easily read and plot
FLurb <- fortify(FLurb)
# class(FLpopdensity)
# str(FLpopdensity)

# plot the shapefile of urban areas onto the ggmap
# what happens when the group argument is removed?
ggmap(google) + geom_polygon(aes(x = long, y = lat, group = group), color='orange', fill = "orange",
                             size=.2, data = FLurb, alpha = 0.8) + 
  geom_point(aes(x = lon, y = lat), color = heat.colors(20, alpha = 0.2)[4], 
             fill = heat.colors(20, alpha = 0.2)[4], size = 1, shape = 21, 
             data = mopa7_FL)

# the ggmap is plotting to the RStudio graphics device
# to save these maps as files in your working directory, you can include ggsave() or tiff()
# above the code that creates the map

# now let's plot maps with and without MOPA sightings side by side
# gridExtra is just one way that worked for my ggmaps, there are many others 
# look up split.screen() and par(mfrow), among others
# facet_wrap is specific to ggplot, in which each level of a given factor can become a different
# facet or panel of the image file, should work with ggmaps
# this arrangement can also be made much prettier with some playing around

library(gridExtra)

plot1 <- ggmap(google, legend = "none") + 
  geom_polygon(aes(x = long, y = lat, group = group), color='orange', fill = "orange",
                                                       size=.2, data = FLurb, alpha = 0.8) + 
  ggtitle("Urban centers in Florida")

plot2 <- ggmap(google, legend = "none") + 
  geom_polygon(aes(x = long, y = lat, group = group), color='orange', fill = "orange",
                                                       size=.2, data = FLurb, alpha = 0.8) + 
  geom_point(aes(x = lon, y = lat), color = heat.colors(20, alpha = 0.2)[4], 
             fill = heat.colors(20, alpha = 0.2)[4], size = 1, shape = 21, 
             data = mopa7_FL) +
  ggtitle("Urban centers and \n MOPA sightings in Florida")

grid.arrange(plot1, plot2, ncol = 2) 


# let's get rid of the city names now, using a different google map type

# download a satellite google map
map <- get_googlemap("florida", zoom = 7, maptype = c("satellite"))

plot1 <- ggmap(map, legend = "none") + 
  geom_polygon(aes(x = long, y = lat, group = group), color='orange', fill = "orange",
                                                       size=.2, data = FLurb, alpha = 0.8) +
  ggtitle("Urban centers in Florida")

plot2 <- ggmap(map, legend = "none") + 
  geom_polygon(aes(x = long, y = lat, group = group), color='orange', fill = "orange",
                                                       size=.2, data = FLurb, alpha = 0.8) + 
  geom_point(aes(x = lon, y = lat), color = heat.colors(20, alpha = 0.2)[4], 
             fill = heat.colors(20, alpha = 0.2)[4], size = 1, shape = 21, 
             data = mopa7_FL) +
  ggtitle("Urban centers and \n MOPA sightings in Florida")

grid.arrange(plot1, plot2, ncol = 2)











