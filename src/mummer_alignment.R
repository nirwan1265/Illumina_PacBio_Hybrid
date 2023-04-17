# Download the latest mummer version
# In vcl command line:
# > wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
# > tar -xzf mummer-4.0.0rc1.tar.gz
# > cd mummer-4.0.0rc1
# > ./configure --prefix=$HOME
# > make isntall
# > make
./nucmer --maxmatch /home/ntanduk/CS_simulated_data/reference/simulated_genome_500k_AT_200_50.fa /home/ntanduk/SPAdes-3.15.5-Linux/bin/results/AT_200_50/Illumina_AT_200_50_pacbio_assembly/scaffolds.fasta -p AT_200_50_alignment

# Alignment of assembled genome to reference genome
./nucmer --maxmatch /Users/nirwan/Desktop/CS_simulated_data/reference/simulated_genome_500k_AT_200_50.fa /Users/nirwan/Desktop/CS_simulated_data/AT_200_50/Illumina_AT_200_50_pacbio_assembly/scaffolds.fasta -p AT_200_50_alignment

# Filter the alignment to keep only 1-to-1 mapping 
./delta-filter -1 AT_200_50_alignment.delta > AT_200_50_alignment_1to1.delta

#make graph
./mummerplot --large --layout --filter --postscript AT_200_50_alignment_1to1.delta -p AT_200_50_alignment.1to1
# install gnuplot
#brew install gnuplot
gnuplot AT_200_50_alignment.1to1.gp
ps2pdf AT_200_50_alignment.1to1.ps AT_200_50_alignment.1to1.pdf

open the pdf file