# R Test Script
# Ada Madejska
# March 1st, 2018
# The purpose of this script is tro learn how to work in R

count <- 0
primes <- c(1,3,7,11)
Names <-
organism <- c("Human", "Chimpanzee", "Yeast")
chromosomes <- c(23, 24, 16)
multicellular <- c(TRUE, TRUE, FALSE)

OrganismTable <- data.frame(organism, chromosomes, multicellular)

write.table(OrganismTable, file = "MyData.csv", row.names = FALSE, na="", col.names = FALSE, sep="")
