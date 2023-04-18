generate_random_genome_sequence <- function(l, r, m, n) {
  # Create a vector of all possible DNA bases
  all_bases <- c("A", "C", "G", "T")
  
  # Generate a vector of random DNA bases
  genome_sequence <- sample(all_bases, l, replace = TRUE)
  
  # Determine the indices to insert the repeat
  repeat_indices <- sample(1:(l - (n * length(r) * m)), n, replace = FALSE)
  repeat_indices <- sort(repeat_indices)
  
  # Insert the repeat at the specified indices
  for(i in 1:n) {
    repeat_start_index <- repeat_indices[i] + ((i - 1) * length(r) * m)
    repeat_end_index <- repeat_start_index + (length(r) * m) - 1
    genome_sequence[repeat_start_index:repeat_end_index] <- rep(r, each = m)
  }
  
  # Collapse the vector into a single string
  genome_sequence_string <- paste(genome_sequence, collapse = "")
  
  return(genome_sequence_string)
}

random_seq <- generate_random_genome_sequence(200000, "AT", 10000, 5)


setwd("~/Desktop/simulated_genome/reference")
# Sequence name and description
seqname <- "Simulated_genome"

library(seqinr)
# Saving a fasta file
write.fasta(sequences = random_seq, names = seqname, file.out = "simulated_genome_200k_AT_10000_5.fa")


generate_random_genome_sequence(l, r, m, n)
l = (numeric) length of the genome, e.g. 500000
r = (character) repeat string, e.g. "AT" or "ATTGC"
m = (numeric) length of the repeat, e.g. 3, so if you do r = "AT", and m = 3 your repeat will be "ATATAT"
n = (numeric) times you will see the repeats, e.g. n = 4 you will se "ATATAT" 4 times randomly in your genome
