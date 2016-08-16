# November  2015
##G.Hurtado
## Basics of R, BGSO NMSU
####################################################################
# Objectives for this portion of the workshop

###Data exploration ######

# 1) histograms
# 2) qq-plots
# 3) scatter plots
# 4) boxplots
# 4) t-tests
# 5) Wilcoxon/ Mann-Whitney

####################################################################


##########################################################3


# upload example data set (BGSO.R.Data) and open the example data set in Excel

setwd("~/Desktop/BGSO_R_workshop_8nov15/")
BGSO.R.Data <- read.csv("BGSO.R.Data.csv", header = TRUE, row.names = NULL)


file.choose() 
# or import data set from environment

# navigate to where the data file is saved (or I you can get it from the thrumb drive or I can send it via email)

# Now we want to get an idea of what the data looks like
# we will use some commands to look at the data
# the commands are: head, tail, and summary
# look at the help for these commands

head (BGSO.R.Data)

tail (BGSO.R.Data)

names (BGSO.R.Data)

# lets start with a histogram of Mass
# the code for a basic histogram is

hist (BGSO.R.Data$mass)

# lets clean this histogram up

# lets start with the title with the term main
?hist
hist (BGSO.R.Data$mass, main = "Mass of K-rats")
# you should see the tile change

# now lets add labels with ylab and xlab terms
hist (BGSO.R.Data$mass, main = "Mass of K-rats", ylab = "Frequency of Mass", xlab ="Mass of K-rats")

# lets change x to 0 to 60  with the xlim term
hist (BGSO.R.Data$mass, main = "Mass of K-rats", ylab = "Frequency of Mass", xlab ="Mass of K-rats", xlim=c(0, 60))

#Lets do another histogram
#lets create a histogram for pter.infect
#lets use the code similar to the final version that we used 
#for the last histogram but change names and labels and change y axis and x axis

hist (BGSO.R.Data$pter.infect, main = "Pter Infections per Animal", ylab = "Frequency of Pter Infections", xlab ="Infections in K-rats", xlim=c(0, 10), ylim=c(0, 100))

###################################
# lets do a normality test for Mass

shapiro.test(BGSO.R.Data$mass)


#Lets make a qq plot from the mass and foot data

qqplot(BGSO.R.Data$mass, BGSO.R.Data$foot)

# Go ahead and modify the qqplot with a title and y and x labels and x, y limits

qqplot(BGSO.R.Data$mass, BGSO.R.Data$foot, main = "qq-plot mass and foot length", ylab = "Foot length", xlab ="Mass of K-rats", xlim=c(0, 60), ylim=c(0, 40))



##################################
# Lets create a scatter plot of mass and foot

plot(BGSO.R.Data$mass, BGSO.R.Data$pter.infect)

# # Go ahead and modify the qqplot with a title and y and x labels and x, y limits

# Now lets modify other parts of the plot- change size of the circle cex term

plot(BGSO.R.Data$mass, BGSO.R.Data$pter.infect, main = "scatter plot", ylab = "Foot length", xlab ="Mass of K-rats", xlim=c(0, 60), ylim=c(0, 10), cex= 2)

# Now change size of title with cex.main term
plot(BGSO.R.Data$mass, BGSO.R.Data$pter.infect, main = "scatter plot", ylab = "Foot length", xlab ="Mass of K-rats", xlim=c(0, 60), ylim=c(0, 10), cex= 2, cex.main=2)

# Change labels and y, x
plot(BGSO.R.Data$mass, BGSO.R.Data$pter.infect, main = "scatter plot", ylab = "Foot length", xlab ="Mass of K-rats", xlim=c(0, 60), ylim=c(0, 10), cex= 2, cex.main=2, cex.lab =2,cex.axis= 1.5)


plot(BGSO.R.Data$mass, BGSO.R.Data$pter.infect, main = "scatter plot", ylab = "Foot length", xlab ="Mass of K-rats", xlim=c(0, 60), ylim=c(0, 10), cex= 2, cex.main=2, cex.lab =2,cex.axis= 1.5)
# can change fonts the same way with font agument e.g. font.main, font.lab...
# can change color with col term and similar agruments as cex and font....
# can use pch...

###############################################

#Lets creat a boxplot to look at mass

boxplot(BGSO.R.Data$mass)

#Now lets do a t-test on mass
t.test(BGSO.R.Data$mass)

t.test(BGSO.R.Data$mass, mu= 0, alt="less", conf=0.95)

#######################################################
#Lets take alook at mass by habitat with boxplots
boxplot(BGSO.R.Data$mass)

boxplot(BGSO.R.Data$mass~ BGSO.R.Data$habitat)

# go ahead and change the boxplot add labels, title .......
### Insert your code here....



###################################
#Wilcox test

#mu= location shift not equal to zero
#alternative to mu = two sided
# nonparametric confidence interval returned in R, conf.int = T
# level of confidence, conf.level = 0.95
# not paried data, paired = F
# to get exact p-value computed, exact = T

wilcox.test(BGSO.R.Data$mass~ BGSO.R.Data$habitat, mu= 0, alt= "two.sided", conf.int= T, conf.level= 0.95, paired= F, exact= T, correct= T)

###########################################
#*Hopefully you found this R workshop helpful!*

