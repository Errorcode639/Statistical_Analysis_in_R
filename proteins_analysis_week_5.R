# R Protein Structure, Expression, and Interactions Script
# Ada Madejska
# March 4th, 2018
# The purpose of this sctipt is to learn about Mass Spec analysis in R

library(ggplot2)
library(gplots)
library(MASS)
library(timeDate)
library(timeSeries)


# Get the data and look at it
peptides.txt <- read.table("peptidefrags.txt", header=FALSE)
peptides <- as.vector(peptides.txt$V1)

# Let's see if the data has a normal distribution
hist(peptides, breaks=400)

mascot.txt <- read.table("mascot.txt", header=FALSE)
mascot <- as.vector(mascot.txt$V1)
xtandem.txt <- read.table("xtandem.txt", header=FALSE)
xtandem <- as.vector(xtandem.txt$V1)
protpro.txt <- read.table("protpro.txt", header=FALSE)
protpro <- as.vector(protpro.txt$V1)

# Merge the data
combinedMSsata <- list(Mascot=mascot, XTandem=xtandem, ProtPro = protpro)
venn(combinedMSsata)

# Load data and limit it to the actual numeric data columns
Dataset <- read.csv("ms.csv", header=TRUE, na.strings="NA", dec=".", strip.white=TRUE)
RawData <- Dataset[,2:14]
filledcols = colSds(RawData) != 0.0
RawData <- RawData[,filledcols]

# Create LDA probabilities for each class X1 based on the values in Dataset
test1.lda <- lda(Dataset$X1~., data=Dataset)
test1.lda.values <- predict(test1.lda, Dataset)

# Take specific columns and use the values as x, y coordinates for plotting
x <- test1.lda.values$x[,1]
y <- test1.lda.values$x[,2]

class <- Dataset$X1
plotData <- data.frame(class, x, y)
centroids <- aggregate(cbind(x,y)~class, plotData, mean)
# Generate distance matrix  from the centroids dataframe
CentrioidDistances <- dist(centroids, method="euclidean", diag=TRUE, upper=FALSE, p=2)

attr(CentrioidDistances, "Labels") <- centroids$class

plot1 <- ggplot(plotData,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids, size=7) + 
  geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") +
  geom_text(aes(label=Dataset$X1), hjust=0, vjust=0, colour="black")
plot1 
