# Simulate a genome
set.seed(123)

# Define repeat sequences
random_nucleotides <- c("A", "T", "C", "G")
direpeats <- paste0(rep(c("AT", "TC", "CG", "GA"), each = 1), collapse = "")
trirepeats <- paste0(rep(c("ATC", "CGT", "TAC", "GCA"), each = 1), collapse = "")
quadrepeats <- paste0(rep(c("ATCG", "CGTA", "GTAC", "TGCA"), each = 1), collapse = "")
longdirepeat <- paste0(rep(c("GT"), times = 60), collapse = "")
#nonrepeat <- "N"

# Define frequencies of each repeat type
#random_nucleotide_freq <- 0.4
direpeat_freq <- 0.1
trirepeat_freq <- 0.1
quadrepeat_freq <- 0.1
longdirepeat_freq <- 0.01
random_nucleotide_freq <- 1 - random_nucleotide_freq - direpeat_freq - trirepeat_freq - quadrepeat_freq - longdirepeat_freq

# Simulating the genome based on the frequencies
repeats <- c(random_nucleotides, direpeats, trirepeats, quadrepeats, longdirepeat)
repeat_freqs <- c(rep(random_nucleotide_freq/4, 4), direpeat_freq, trirepeat_freq, quadrepeat_freq, longdirepeat_freq)
seq <- sample(repeats, 500000, replace = TRUE, prob = repeat_freqs)
seq <- paste(seq, collapse = "")
seq




# Sequence name and description
seqname <- "Simulated_genome"

# Saving a fasta file
write.fasta(sequences = seq, names = seqname, file.out = "simulated_reference/simulated_genome.fa")

