# R Genomic analysis Test Script
# Ada Madejska
# March 3rd, 2018
# The purpose of this script is to find ORFs using R

library(Biostrings)
library(seqinr)

# Get the seq
AB003468 <- readDNAStringSet("AB003468.fasta")
AB003468 = as.character(AB003468)
sequence <- AB003468

# Create start and stop variables of DNA
start_codon <- "ATG"
stop_codons <- c("TAG", "TAA", "TGA")

# Create variables to store the results in
start_pos <- c()
revstart_pos <- c()
stop_pos <- c()
revstop_pos <- c()

matches <- matchPattern(start_codon, sequence)

# Take staring location of each match and store it
start_pos <- c(start_pos, start(matches))

# Do the same for the reverse string
revmatches <- matchPattern(reverseComplement(DNAString(start_codon)), sequence)
revstart_pos <- c(revstart_pos, start(revmatches))

# Sort outpouts numerically
start_pos <- sort(start_pos)
revstart_pos <- sort(revstart_pos, decreasing=TRUE)

# Find STOP codons
for (codon in stop_codons) {
  matches <- matchPattern(codon, sequence)
  stop_pos <- c(stop_pos, start(matches))
  revmatches <- matchPattern(reverseComplement(DNAString(codon)), sequence)
  revstop_pos <- c(revstop_pos, start(revmatches))
}

stop_pos <- sort(stop_pos)
revstop_pos <- sort(revstart_pos, decreasing=TRUE)

# Now we need to map out the actual ORFs that might exist in the sequence
k <- 150 # Threshold size of ORF
lengths <- vector(mode="numeric")
stop_pointers <- c(0,0,0) # Will hold the location of each stop in each readingframe
count <- 0

for (current_start in start_pos) {
  frame <- (current_start%%3) + 1
  stop_pointer <- stop_pointers[frame]
  if (stop_pointer <= length(stop_pos) && (stop_pointer == 0 || stop_pos[stop_pointer] < current_start)) {
    stop_pointer <- stop_pointer + 1
    while ((stop_pointer<= length(stop_pos)) && ((stop_pos[stop_pointer] <= current_start) 
                                                 || (((stop_pos[stop_pointer]%%3)+ 1) != frame ))) {
      stop_pointer <- stop_pointer + 1
    }
    stop_pointers[frame] <- stop_pointer
    if (stop_pointer <= length(stop_pos)) {
      if((stop_pos[stop_pointer]+2 - current_start+1) > k) {
        count <- count + 1
        print(count)
        print("Frame:")
        print(frame)
        print("Start:")
        print(current_start)
        print("Stop:")
        print(stop_pos[stop_pointer])
        print("Length:")
        print(stop_pos[stop_pointer]+2-current_start+1)
        lengths <- c(lengths, (stop_pos[stop_pointer]+2-current_start+1))
        print("Sequence:")
        print(subseq(sequence, current_start, stop_pos[stop_pointer]+2 ))
      }
    }
  }
}

lengths <- sort(lengths)

# Create a histogram to visualize the results
bins <- seq(0,1000,50)
hist(lengths, breaks=bins, col="red", xlim=c(0,1000))