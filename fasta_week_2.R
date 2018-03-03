# R FASTA Test Script
# Ada Madejska
# March 2nd, 2018
# The puropse of this script is to learn about using FASTA files in R

library(ape)
library(seqinr)

# Read a local FASTA file
cox1 <- read.fasta(file="cox1.fasta", seqtype="AA")
seq1 <- cox1[1]

# Download a cloning vector from GenBank as binary and save it as FASTA
AB003468 <- read.GenBank("AB003468", as.character="TRUE")
write.dna(AB003468, file="AB003468.fasta", format="fasta", append=FALSE, nbcol=6, colsep=" ", colw=10)

# Get the sequence of Ab003468 and count the number the nucleotides
CloningVector <- AB003468[[1]]
count <- count(CloningVector,1)
count

# Inspect the ratio of G/C residues compared to A/T
GC <- GC(CloningVector)
GC

# Break the sequence into 200 long chunks and calculate GC content for each section
GCwindow <- seq(1, length(CloningVector)-200, by=200)
n <- length(GCwindow)
Chunks <- numeric(n)
for (i in 1:n) {
  chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  chunkGC <- GC(chunk)
  print(chunkGC)
  Chunks[i] <- chunkGC
}

# Plot the output to see how the GC content changes over the sequence
plot(GCwindow, Chunks, type="b", xlab="Nucleotide start position", ylab="GC content")

# Create a function for customizable GC analysis
slidingwindowGCplot <- function(windowsize, inputseq) {
  GCwindow <- seq(1, length(inputseq)-windowsize, by=windowsize)
  
  # Find the length of GCwindow
  n <- length(GCwindow)
  
  # Make a vector of the same size as GCwindow length
  Chunks <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    Chunks[i] <- chunkGC
  }
  
  # Plot the output to see how the GC content changes over the sequence
  plot(GCwindow, Chunks, type="b", xlab="Nucleotide start position", ylab="GC content", main=paste("GC Plot with windowsize ", windowsize))
}
