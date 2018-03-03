# R Test Script
# Ada Madejska
# March 1st, 2018
# The purpose of this script is tro learn how to work in R

primes <- c(1,3,7,11)
Names <- c("Bob", "Ted", "Carol", "Alice")
Truth <- c(TRUE, FALSE)

organism <- c("Human", "Chimpanzee", "Yeast")
chromosomes <- c(23, 24, 16)
multicellular <- c(TRUE, TRUE, FALSE)

OrganismTable <- data.frame(organism, chromosomes, multicellular)

write.table(OrganismTable, file = "MyData.csv", row.names = FALSE, na="", col.names = TRUE, sep="")

NewDataFrame <- read.csv("MyData.csv", header=TRUE, sep="")

# Check how many organisms have more than 20 chromosomes
count <- 0
for (val in OrganismTable$chromosomes) {
  if(val>20)  count = count + 1
}
print(count)

barplot(OrganismTable$chromosomes)
