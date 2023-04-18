# Running spades
# git clone

#!/bin/bash

# set the input and output directory paths
illumina_dir="illumina"
pacbio_dir="pacbio"
reference_dir="reference"
output_dir="spades_assemblies"

# loop through each reference file in the input directory
for ref_file in "${reference_dir}"/*.fa
do
# get the filename without the path or extension
filename=$(basename "${ref_file}")
filename="${filename%.*}"

# run the command with the current input files and output directory
# -t = number of threads. in vcl it is 16
spades.py --isolate -t 16 \
-1 "${illumina_dir}/illumina_simulated_l50_c30_${filename}_pe1.fq" \
-2 "${illumina_dir}/illumina_simulated_l50_c30_${filename}_pe2.fq" \
--pacbio "${pacbio_dir}/simulated_pacbio_c20_l12000_${filename}__0001.fastq" \
--trusted-contigs "${reference_dir}/${filename}.fa" \
-o "${output_dir}/${filename}"
done


#./spades.py --isolate -t 16 -1 /home/ntanduk/CS_simulated_data/illumina/illumina_simulated_l50_c30_simulated_genome_500k_AT_1500_50_pe1.fq -2 /home/ntanduk/CS_simulated_data/illumina/illumina_simulated_l50_c30_simulated_genome_500k_AT_1500_50_pe2.fq --pacbio /home/ntanduk/CS_simulated_data/pacbio/simulated_pacbio_c20_l12000_simulated_genome_500k_AT_1500_50__0001.fastq --trusted-contigs /home/ntanduk/CS_simulated_data/reference/simulated_genome_500k_AT_1500_50.fa -o results/Illumina_AT_1500_50_pacbio_assembly

