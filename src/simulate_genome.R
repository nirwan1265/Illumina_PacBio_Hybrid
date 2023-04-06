# Simulate a genome
set.seed(123)

# Define repeat sequences
monorepeats <- paste0(rep(c("A", "T", "C", "G"), each = 1), collapse = "")
direpeats <- paste0(rep(c("AT", "TC", "CG", "GA"), each = 1), collapse = "")
trirepeats <- paste0(rep(c("ATC", "CGT", "TAC", "GCA"), each = 1), collapse = "")
drepeat <- paste0(rep(c("ATCG", "CGTA", "GTAC", "TGCA"), each = 1), collapse = "")
longrepeat <- paste0(rep(c("GGGGGGGGGG", "AAAAAAAAAAAAAAAA","CCCCCCCCCCC","TTTTTTTTTTT"), each = 1), collapse = "")

# Define frequencies of each repeat type
monorepeat_freq <- 0.01
direpeat_freq <- 0.01
trirepeat_freq <- 0.01
drepeat_freq <- 0.01
longrepeat_freq <- 0.01
nonrepeat_freq <- 1 - monorepeat_freq - direpeat_freq - trirepeat_freq - drepeat_freq - longrepeat_freq

# Simulating the genome based on the frequencies
repeats <- c(monorepeats, direpeats, trirepeats, drepeat, longrepeat)
repeat_freqs <- c(monorepeat_freq, direpeat_freq, trirepeat_freq, drepeat_freq, longrepeat_freq)
seq <- sample(repeats, 100, replace = TRUE, prob = repeat_freqs)
seq <- paste(seq, collapse = "")
seq <- paste(sample(c(rep("N", 60)), 1), seq, sample(c(rep("N", 60)), 1), sep = "")

# Sequence name and description
seqname <- "Simulated1"

# Saving a fasta file
write.fasta(sequences = seq, names = seqname, file.out = "simulated1.fa")

