# R Sequence Alignment Test Script
# Ada Madejska
# March 3rd, 2018
# The purpose of this script is to learn how to use sequence alignment in R

library(Biostrings)
library(phangorn)
library(seqinr)
library(msa)

# Read the fasta file
prokaryotes <- read.fasta(file="prok.fasta", seqtype="DNA")

# Parse and convert to a simple text first 2 sequences
seq1 <- as.character(prokaryotes[[1]])
seq1 = paste(seq1, collapse="")

seq2 <- as.character(prokaryotes[[2]])
seq2 = paste(seq2, collapse="")

# Align seq1 and seq2 using default settings of Biostrings
pairalign <- pairwiseAlignment(pattern=seq2, subject=seq1)
summary(pairalign)

# Export alignment back to FASTA file
pairalignString = BStringSet(c(toString(subject(pairalign)), toString(pattern(pairalign))))
writeXStringSet(pairalignString, "aligned.txt", format="FASTA")

# Load cox data and put each sequence in separate variable
coxgenes <- read.fasta(file="cox1multi.fasta", seqtype="AA")
cox1 <- as.character(coxgenes[[1]])
cox2 <- as.character(coxgenes[[2]])

# Generate simple dot plot of cox data
dotPlot(cox1, cox2, wsize=3, wstep=3, nmatch=3, main="Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Read fasta files as character sets
coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

# Align cox1 protein sequence using CLUSTALW ( used to look at the similarities/differences btwn groups of sequences)
coxAligned <- msa(coxAA)

prokAligned <- msa(prokDNA)
print(prokAligned, show="complete")

#Export msa 
prokAlignStr = as(prokAligned, "DNAStringSet")
writeXStringSet(prokAlignStr, file="prokAligned.fasta")

coxAlignStr = as(coxAligned, "AAStringSet")
writeXStringSet(coxAlignStr, file="coxAligned.fasta")

# Write into a PHILIP Format
write.phylip(coxAligned, "coxAligned.phylip")

# Convert prokAligned to seqinr format
prokAligned2 <- msaConvert(prokAligned, type="seqinr::alignment")

# Generate a distance matrix
prokdist <- dist.alignment(prokAligned2, "identity")