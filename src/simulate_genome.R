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

random_seq <- generate_random_genome_sequence(500000, "ATTGGTTA", 1000, 200)


setwd("~/Desktop/simulated_genome/reference")
# Sequence name and description
seqname <- "Simulated_genome"

# Saving a fasta file
write.fasta(sequences = random_seq, names = seqname, file.out = "simulated_genome_500k_ATTGGTTA_1000_200.fa")

