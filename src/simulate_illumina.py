
# First create conda environment
# > conda create -n reseq_simulate_illumina
# > conda activate reseq_simulate_illumina

# Second install the reseq
# > conda install -c bioconda -c conda-forge reseq

# installation protocols followed using this website:
# https://github.com/scchess/Art/blob/master/art_illumina_README

# 1st step install the binary version of the software (you can compile it yourself too)
# https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

# 2nd step: add GSL libraries (in command line)
# > export CFLAGS="$CFLAGS -I/opt/local/include" CPPFLAGS="$CPPFLAGS -I/opt/local/include" LDFLAGS="$LDFLAGS -L/opt/local/lib"
# > export CFLAGS="$CFLAGS -I/usr/local/include" CPPFLAGS="$CPPFLAGS -I/usr/local/include" LDFLAGS="$LDFLAGS -L/usr/local/lib"

# 3rd step: compile and install
# go to that folder where you downloaded the binary version of the software
# > ./configure --prefix=$HOME
# > make 
# > make install
# ./art_illumina

# to run the simulation
./art_illumina -ss NS50 -sam -i simluated_genome/simulated_genome_100kb.fa -l 50 -m 60 -s 10 -f 10 -p -o illumina_simulated_500kb_l50_c10_pe

# -ss is the type of illumina machine
# i is the simulated reference file
# -l is the length of the reads
# -f is the coverage
# - o is the output
# -p for paired ends
# -m the mean size of DNA/RNA fragments for paired-end simulations

