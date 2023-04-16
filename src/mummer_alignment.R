# Download the latest mummer version
# In vcl command line:
# > wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
# > tar -xzf mummer-4.0.0rc1.tar.gz
# > cd mummer-4.0.0rc1
# > ./configure --prefix=$HOME
# > make isntall
# > make
nucmer --maxmatch /home/ntanduk/CS_simulated_data/reference/simulated_genome_500k_AT_200_50.fa /home/ntanduk/SPAdes-3.15.5-Linux/bin/results/AT_200_50/Illumina_AT_200_50_pacbio_assembly/scaffolds.fasta -p AT_200_50_alignment