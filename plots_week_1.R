# R Test Script
# Ada Madejska
# March 2nd, 2018
# The purpose of this script is tro learn how to work in R: plots

library(ggplot2)
library(reshape2)

# Obtain data from a csv file
rawdata <- read.csv("Week_1_Plotdata.csv", header=TRUE)
ggplot(rawdata, aes(x=Subject, y=a)) + geom_point() # plot rawdata using ggplot

# Convert  the data to a simple X Y coordinates by Subject  and condition
melted = melt(rawdata, id.vars="Subject", measure.vars=c("a", "c", "d", "e", "f", "g", "j", "k"))

# Create a plot from the melted data frame
myPlot <- ggplot(melted, aes(x=variable, y=value, col=Subject, group=Subject)) + geom_point() + geom_line() +
  xlab("Sample") +
  ylab("# Observed") + 
  ggtitle("Some observations I made in the lab")

# Display myPlot
myPlot

# Save myPlot to pdf
ggsave(filename="MyPlot.pdf", plot=myPlot)