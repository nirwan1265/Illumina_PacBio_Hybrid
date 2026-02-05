# Download quast:
# on Mac or linux, for windows check this: https://quast.sourceforge.net/docs/manual.html#sec1
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
tar -xzf quast-5.2.0.tar.gz
cd quast-5.2.0

# run genome assembly comparison 
quast.py -t 10 -r reference_genome.fasta -o output_dir contigs.fasta scaffolds.fasta



#
./quast.py -t 1 -r simulated_genome_10k_AT_1000_5.fa -o genome_10k_AT_1000_5_illuminaonly contigs_10k_AT_1000_5_illuminaonly.fasta scaffolds_10k_AT_1000_5_illuminaonly.fasta
